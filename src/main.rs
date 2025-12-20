// Copyright (c) 2024 Maxwell Campbell. Licensed under the MIT License.
use clap::error::ErrorKind;
use clap::{CommandFactory, Parser};
use indicatif::{ProgressBar, ProgressStyle};
use pdbtbx::{PDBError, ReadOptions};
use quick_xml::SeError as XmlError;
use rust_sasa::options::SASACalcError;
use rust_sasa::options::SASAOptions;
use rust_sasa::structures::atomic::SASAResult;
use rust_sasa::utils::configure_thread_pool;
use rust_sasa::{sasa_result_to_json, sasa_result_to_protein_object, sasa_result_to_xml};
use snafu::{ResultExt, Snafu};

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

#[derive(clap::ValueEnum, Clone, Default, Debug)]
pub enum SASALevel {
    Atom,
    #[default]
    Residue,
    Chain,
    Protein,
}

#[derive(clap::ValueEnum, Clone, Default, Debug)]
enum OutputFormat {
    Xml,
    #[default]
    Json,
    Pdb,
    Cif,
}

impl OutputFormat {
    fn file_extension(&self) -> &'static str {
        match self {
            OutputFormat::Xml => "xml",
            OutputFormat::Json => "json",
            OutputFormat::Pdb => "pdb",
            OutputFormat::Cif => "cif",
        }
    }

    fn from_file_extension(filename: &str) -> Self {
        match filename.rsplit('.').next() {
            Some("xml") => OutputFormat::Xml,
            Some("json") => OutputFormat::Json,
            Some("pdb") => OutputFormat::Pdb,
            Some("cif") => OutputFormat::Cif,
            _ => OutputFormat::Json, // default
        }
    }
}

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// File to read from.
    #[arg()]
    input: String,

    /// Output file path.
    #[arg()]
    output: String,

    /// Output depth. (i.e: protein, chain, residue, atom)
    #[arg(short, long, default_value_t, value_enum)]
    output_depth: SASALevel,

    /// Output format (required when processing directories, inferred from file extension for single files)
    #[arg(short, long, value_enum)]
    format: Option<OutputFormat>,

    /// Number of Shrake Rupley points
    #[arg(short, long, default_value_t = 100)]
    n_points: usize,

    /// Probe radius in Angstroms (default: 1.4)
    #[arg(short, long, default_value_t = 1.4)]
    probe_radius: f32,

    /// Include hydrogen atoms in SASA calculation (default: hydrogens are excluded)
    #[arg(short = 'H', long, default_value_t = false)]
    include_hydrogens: bool,

    /// Path to custom radii configuration file (default: uses embedded protor.config)
    #[arg(short = 'r', long)]
    radii_file: Option<String>,

    /// Allow fallback to van der Waals radii when custom radius is not found (default: false, strict mode)
    #[arg(short = 'a', long, default_value_t = false)]
    allow_vdw_fallback: bool,

    /// Include HETATM records. Defaults to false as ProtOr config does not contain records for non-standard amino acids. (default: false)
    #[arg(short = 'n', long, default_value_t = false)]
    include_hetatms: bool,

    /// Configure number of threads used to parallelize SASA computation. Mode: -1 uses all CPU cores. (default: -1)
    #[arg(short = 't', long, default_value_t = -1)]
    threads: isize,

    #[arg(short = 'R', long, default_value_t = false)]
    read_radii_from_occupancy: bool,
}

#[derive(Debug, Snafu)]
pub enum CLIError {
    #[snafu(display("SASA calculation failed: {source}"))]
    SASACalculation { source: SASACalcError },

    #[snafu(display("Failed to create thread pool: {source}"))]
    ThreadPool { source: std::io::Error },

    #[snafu(display("Failed to read from input file"))]
    InputFileRead { errors: Vec<PDBError> },

    #[snafu(display("Failed to serialize to XML: {source}"))]
    XMLSerialization { source: XmlError },

    #[snafu(display("Failed to serialize to JSON: {source}"))]
    JSONSerialization { source: serde_json::Error },

    #[snafu(display("Failed to serialize to Protein Object: {message}"))]
    ProteinSerialization { message: String },

    #[snafu(display("Failed to write PDB"))]
    PDBWrite { errors: Vec<PDBError> },

    #[snafu(display("Failed to write MMCIF"))]
    MMCIFWrite { errors: Vec<PDBError> },

    #[snafu(display("Failed to read directory: {source}"))]
    DirectoryRead { source: std::io::Error },

    #[snafu(display("Failed to write output file: {source}"))]
    FileWrite { source: std::io::Error },

    #[snafu(display("Failed to load radii file: {source}"))]
    RadiiFileLoad { source: std::io::Error },

    #[snafu(display("Input path does not exist: {path}"))]
    InputPathNotFound { path: String },

    #[snafu(display("Input path appears to be a directory but does not exist: {path}"))]
    InputDirectoryNotFound { path: String },
}

impl CLIError {
    fn to_clap_error(&self) -> clap::Error {
        let msg = match self {
            Self::InputFileRead { errors }
            | Self::PDBWrite { errors }
            | Self::MMCIFWrite { errors } => {
                let error_details = errors
                    .iter()
                    .map(|e| e.to_string())
                    .collect::<Vec<_>>()
                    .join("; ");
                format!("{self}: {error_details}")
            }
            _ => self.to_string(),
        };

        Args::command().error(ErrorKind::Format, msg)
    }
}

#[allow(clippy::too_many_arguments)]
fn process(
    input_file: String,
    output_file: String,
    output_depth: SASALevel,
    format: &OutputFormat,
    n_points: usize,
    probe_radius: f32,
    threads: isize,
    include_hydrogens: bool,
    radii_file: Option<&str>,
    allow_vdw_fallback: bool,
    include_hetatms: bool,
    read_radii_from_occupancy: bool,
) -> Result<(), CLIError> {
    let (pdb, _) = ReadOptions::default()
        .set_level(pdbtbx::StrictnessLevel::Loose)
        .read(input_file)
        .map_err(|errors| CLIError::InputFileRead { errors })?;
    let result = calculate_sasa_and_wrap(
        &pdb,
        &output_depth,
        n_points,
        probe_radius,
        threads,
        include_hydrogens,
        radii_file,
        allow_vdw_fallback,
        include_hetatms,
        read_radii_from_occupancy,
    )
    .context(SASACalculationSnafu)?;

    match format {
        OutputFormat::Xml => {
            let r = sasa_result_to_xml(&result).context(XMLSerializationSnafu)?;
            std::fs::write(output_file, r).context(FileWriteSnafu)?;
        }
        OutputFormat::Json => {
            let r = sasa_result_to_json(&result).context(JSONSerializationSnafu)?;
            std::fs::write(output_file, r).context(FileWriteSnafu)?;
        }
        OutputFormat::Pdb => {
            let mut original_pdb = pdb.clone();
            sasa_result_to_protein_object(&mut original_pdb, &result)
                .map_err(|message| CLIError::ProteinSerialization { message })?;
            pdbtbx::save(&original_pdb, output_file, pdbtbx::StrictnessLevel::Loose)
                .map_err(|errors| CLIError::PDBWrite { errors })?;
        }
        OutputFormat::Cif => {
            let mut original_pdb = pdb.clone();
            sasa_result_to_protein_object(&mut original_pdb, &result)
                .map_err(|message| CLIError::ProteinSerialization { message })?;
            pdbtbx::save(&original_pdb, output_file, pdbtbx::StrictnessLevel::Loose)
                .map_err(|errors| CLIError::MMCIFWrite { errors })?;
        }
    }
    Ok(())
}

/// Macro to reduce duplication in calculate_sasa_and_wrap
macro_rules! process_level {
    ($level_constructor:ident, $result_variant:ident, $pdb:expr, $n_points:expr, $probe_radius:expr, $threads:expr, $include_hydrogens:expr, $radii_file:expr, $allow_vdw_fallback:expr, $include_hetatms:expr, $read_radii_from_occupancy:expr) => {{
        let mut options = SASAOptions::$level_constructor()
            .with_threads($threads)
            .with_n_points($n_points)
            .with_probe_radius($probe_radius)
            .with_include_hydrogens($include_hydrogens)
            .with_allow_vdw_fallback($allow_vdw_fallback)
            .with_include_hetatms($include_hetatms)
            .with_read_radii_from_occupancy($read_radii_from_occupancy);
        if let Some(path) = $radii_file {
            options = options
                .with_radii_file(path)
                .map_err(|e| SASACalcError::RadiiFileLoad { source: e })?;
        }
        let result = options.process($pdb)?;
        Ok(SASAResult::$result_variant(result))
    }};
}

#[allow(clippy::too_many_arguments)]
fn calculate_sasa_and_wrap(
    pdb: &pdbtbx::PDB,
    level: &SASALevel,
    n_points: usize,
    probe_radius: f32,
    threads: isize,
    include_hydrogens: bool,
    radii_file: Option<&str>,
    allow_vdw_fallback: bool,
    include_hetatms: bool,
    read_radii_from_occupancy: bool,
) -> Result<SASAResult, SASACalcError> {
    match level {
        SASALevel::Atom => process_level!(
            atom_level,
            Atom,
            pdb,
            n_points,
            probe_radius,
            threads,
            include_hydrogens,
            radii_file,
            allow_vdw_fallback,
            include_hetatms,
            read_radii_from_occupancy
        ),
        SASALevel::Residue => process_level!(
            residue_level,
            Residue,
            pdb,
            n_points,
            probe_radius,
            threads,
            include_hydrogens,
            radii_file,
            allow_vdw_fallback,
            include_hetatms,
            read_radii_from_occupancy
        ),
        SASALevel::Chain => process_level!(
            chain_level,
            Chain,
            pdb,
            n_points,
            probe_radius,
            threads,
            include_hydrogens,
            radii_file,
            allow_vdw_fallback,
            include_hetatms,
            read_radii_from_occupancy
        ),
        SASALevel::Protein => process_level!(
            protein_level,
            Protein,
            pdb,
            n_points,
            probe_radius,
            threads,
            include_hydrogens,
            radii_file,
            allow_vdw_fallback,
            include_hetatms,
            read_radii_from_occupancy
        ),
    }
}

fn ensure_output_directory(output_path: &str) -> Result<(), CLIError> {
    let output_dir = std::path::Path::new(output_path);

    // If path exists but is not a directory, return an error
    if output_dir.exists() && !output_dir.is_dir() {
        return Err(CLIError::DirectoryRead {
            source: std::io::Error::new(
                std::io::ErrorKind::AlreadyExists,
                "Output path exists but is not a directory",
            ),
        });
    }

    // Create directory if it doesn't exist
    if !output_dir.exists() {
        std::fs::create_dir_all(output_dir).map_err(|e| CLIError::DirectoryRead { source: e })?;
    }

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn process_directory(
    input_dir: &str,
    output_dir: &str,
    output_depth: SASALevel,
    n_points: usize,
    probe_radius: f32,
    format: &OutputFormat,
    include_hydrogens: bool,
    radii_file: Option<&str>,
    allow_vdw_fallback: bool,
    include_hetatms: bool,
    read_radii_from_occupancy: bool,
) -> Result<(), CLIError> {
    use rayon::prelude::*;
    use std::sync::Mutex;

    ensure_output_directory(output_dir)?;

    let errors = Mutex::new(Vec::new());
    let files: Vec<_> = std::fs::read_dir(input_dir)
        .map_err(|e| CLIError::DirectoryRead { source: e })?
        .collect();

    // Create progress bar
    let pb = ProgressBar::new(files.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})",
            )
            .expect("Progress bar template should be valid")
            .progress_chars("#>-"),
    );
    files.par_iter().for_each(|entry| match entry {
        Ok(entry) => {
            let path = entry.path();
            if path.is_file() {
                // Convert path to string, skip if non-UTF8
                let Some(input_path) = path.to_str() else {
                    errors
                        .lock()
                        .expect("Mutex should not be poisoned")
                        .push(format!(
                            "Skipping file with non-UTF8 path: {}",
                            path.display()
                        ));
                    pb.inc(1);
                    return;
                };

                // Extract filename, skip if missing
                let Some(filename_os) = path.file_stem() else {
                    errors
                        .lock()
                        .expect("Mutex should not be poisoned")
                        .push(format!("Skipping file without name: {}", path.display()));
                    pb.inc(1);
                    return;
                };

                let Some(filename) = filename_os.to_str() else {
                    errors
                        .lock()
                        .expect("Mutex should not be poisoned")
                        .push(format!(
                            "Skipping file with non-UTF8 name: {}",
                            path.display()
                        ));
                    pb.inc(1);
                    return;
                };

                let extension = format.file_extension();
                let output_path =
                    std::path::Path::new(output_dir).join(format!("{filename}.{extension}"));

                // Convert output path to string, skip if non-UTF8
                let Some(output_path_str) = output_path.to_str() else {
                    errors
                        .lock()
                        .expect("Mutex should not be poisoned")
                        .push(format!(
                            "Skipping file. Output path is non-UTF8: {}",
                            output_path.display()
                        ));
                    pb.inc(1);
                    return;
                };

                pb.set_message(format!("Processing {filename}"));
                match process(
                    input_path.to_string(),
                    output_path_str.to_string(),
                    output_depth.clone(),
                    format,
                    n_points,
                    probe_radius,
                    1, // Single-threaded for individual files (directory parallelism handled by rayon)
                    include_hydrogens,
                    radii_file,
                    allow_vdw_fallback,
                    include_hetatms,
                    read_radii_from_occupancy,
                ) {
                    Ok(_) => pb.inc(1),
                    Err(e) => {
                        errors
                            .lock()
                            .expect("Mutex should not be poisoned")
                            .push(format!("Error processing {filename}: {e}"));
                        pb.inc(1);
                    }
                }
            }
        }
        Err(e) => {
            errors
                .lock()
                .expect("Mutex should not be poisoned")
                .push(format!("Error reading directory entry: {e}"));
        }
    });

    pb.finish_with_message("Processing complete!");

    // Report errors at the end
    let errors = errors.into_inner().expect("Mutex should not be poisoned");
    if !errors.is_empty() {
        eprintln!("\nThe following errors occurred during processing:");
        for error in &errors {
            eprintln!("  - {error}");
        }
        eprintln!("\nTotal errors: {}", errors.len());
    } else {
        println!("All files processed successfully!");
    }

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn process_single_file(
    input_file: String,
    output_file: String,
    output_depth: SASALevel,
    n_points: usize,
    probe_radius: f32,
    include_hydrogens: bool,
    radii_file: Option<&str>,
    allow_vdw_fallback: bool,
    include_hetatms: bool,
    read_radii_from_occupancy: bool,
    threads: isize,
) -> Result<(), CLIError> {
    println!("Processing single file...");

    // Create parent directory for output file if it doesn't exist
    if let Some(parent) = std::path::Path::new(&output_file).parent() {
        if !parent.as_os_str().is_empty() && !parent.exists() {
            std::fs::create_dir_all(parent).map_err(|e| CLIError::FileWrite { source: e })?;
        }
    }

    let format = OutputFormat::from_file_extension(&output_file);

    process(
        input_file,
        output_file,
        output_depth,
        &format,
        n_points,
        probe_radius,
        threads,
        include_hydrogens,
        radii_file,
        allow_vdw_fallback,
        include_hetatms,
        read_radii_from_occupancy,
    )?;
    println!("Finished!");
    Ok(())
}

fn main() {
    let args = Args::parse();

    if let Err(e) = run(args) {
        e.to_clap_error().exit();
    }
}

fn run(args: Args) -> Result<(), CLIError> {
    let input_path = std::path::Path::new(&args.input);
    let radii_file = args.radii_file.as_deref();

    // Configure thread pool based on user preference
    configure_thread_pool(args.threads).map_err(|e| CLIError::ThreadPool { source: e })?;

    if !input_path.exists() {
        // If path ends with '/' or '\' it is probably meant to be a directory
        if args.input.ends_with('/') || args.input.ends_with('\\') {
            return Err(CLIError::InputDirectoryNotFound {
                path: args.input.clone(),
            });
        }
        // Otherwise throw generic "not found" error
        return Err(CLIError::InputPathNotFound {
            path: args.input.clone(),
        });
    }

    if input_path.is_dir() {
        let format = args.format.ok_or_else(|| CLIError::DirectoryRead {
            source: std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "Format argument is required when processing directories",
            ),
        })?;
        process_directory(
            &args.input,
            &args.output,
            args.output_depth,
            args.n_points,
            args.probe_radius,
            &format,
            args.include_hydrogens,
            radii_file,
            args.allow_vdw_fallback,
            args.include_hetatms,
            args.read_radii_from_occupancy,
        )
    } else {
        process_single_file(
            args.input,
            args.output,
            args.output_depth,
            args.n_points,
            args.probe_radius,
            args.include_hydrogens,
            radii_file,
            args.allow_vdw_fallback,
            args.include_hetatms,
            args.read_radii_from_occupancy,
            args.threads,
        )
    }
}
