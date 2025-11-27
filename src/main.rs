use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use pdbtbx::{PDBError, ReadOptions};
use quick_xml::SeError as XmlError;
use rust_sasa::options::SASACalcError;
use rust_sasa::options::SASAOptions;
use rust_sasa::structures::atomic::SASALevel;
use rust_sasa::structures::atomic::SASAResult;
use rust_sasa::{sasa_result_to_json, sasa_result_to_protein_object, sasa_result_to_xml};
use snafu::{ResultExt, Snafu};

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

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
}

#[derive(Debug, Snafu)]
pub enum CLIError {
    #[snafu(display("Element missing for atom"))]
    SASACalculation { source: SASACalcError },

    #[snafu(display("Failed to read from input file"))]
    InputFileRead { errors: Vec<PDBError> },

    #[snafu(display("Failed to serialize to XML"))]
    XMLSerialization { source: XmlError },

    #[snafu(display("Failed to serialize to JSON"))]
    JSONSerialization { source: serde_json::Error },

    #[snafu(display("Failed to serialize to Protein Object"))]
    ProteinSerialization { message: String },

    #[snafu(display("Failed to write PDB"))]
    PDBWrite { errors: Vec<PDBError> },

    #[snafu(display("Failed to write MMCIF"))]
    MMCIFWrite { errors: Vec<PDBError> },

    #[snafu(display("Failed to read directory"))]
    DirectoryRead { source: std::io::Error },

    #[snafu(display("Failed to write output file"))]
    FileWrite { source: std::io::Error },
}

fn process(
    input_file: String,
    output_file: String,
    output_depth: SASALevel,
    format: &OutputFormat,
    n_points: usize,
    parallel: bool,
) -> Result<(), CLIError> {
    let (pdb, _) = ReadOptions::default()
        .set_level(pdbtbx::StrictnessLevel::Loose)
        .read(input_file)
        .map_err(|errors| CLIError::InputFileRead { errors })?;
    let result = calculate_sasa_and_wrap(&pdb, &output_depth, n_points, parallel)
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

/// Backward compatibility helper function
fn calculate_sasa_and_wrap(
    pdb: &pdbtbx::PDB,
    level: &SASALevel,
    n_points: usize,
    parallel: bool,
) -> Result<SASAResult, SASACalcError> {
    match level {
        SASALevel::Atom => {
            let result = SASAOptions::atom_level()
                .with_parallel(parallel)
                .with_n_points(n_points)
                .process(pdb)?;
            Ok(SASAResult::Atom(result))
        }
        SASALevel::Residue => {
            let result = SASAOptions::residue_level()
                .with_parallel(parallel)
                .with_n_points(n_points)
                .process(pdb)?;
            Ok(SASAResult::Residue(result))
        }
        SASALevel::Chain => {
            let result = SASAOptions::chain_level()
                .with_parallel(parallel)
                .with_n_points(n_points)
                .process(pdb)?;
            Ok(SASAResult::Chain(result))
        }
        SASALevel::Protein => {
            let result = SASAOptions::protein_level()
                .with_parallel(parallel)
                .with_n_points(n_points)
                .process(pdb)?;
            Ok(SASAResult::Protein(result))
        }
    }
}

fn validate_output_directory(output_path: &str) -> Result<(), CLIError> {
    let output_dir = std::path::Path::new(output_path);
    if !output_dir.is_dir() {
        return Err(CLIError::DirectoryRead {
            source: std::io::Error::new(
                std::io::ErrorKind::NotFound,
                "Output path is not a directory",
            ),
        });
    }
    Ok(())
}

fn process_directory(
    input_dir: &str,
    output_dir: &str,
    output_depth: SASALevel,
    n_points: usize,
    format: &OutputFormat,
) -> Result<(), CLIError> {
    use rayon::prelude::*;
    use std::sync::Mutex;

    validate_output_directory(output_dir)?;

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
            .unwrap()
            .progress_chars("#>-"),
    );

    files.par_iter().for_each(|entry| match entry {
        Ok(entry) => {
            let path = entry.path();
            if path.is_file() {
                let input_path = path.to_str().unwrap().to_string();
                let filename = path.file_stem().unwrap().to_str().unwrap();
                let extension = format.file_extension();
                let output_path = format!("{output_dir}/{filename}.{extension}");

                pb.set_message(format!("Processing {filename}"));
                match process(
                    input_path,
                    output_path,
                    output_depth.clone(),
                    format,
                    n_points,
                    false,
                ) {
                    Ok(_) => pb.inc(1),
                    Err(e) => {
                        errors
                            .lock()
                            .unwrap()
                            .push(format!("Error processing {filename}: {e}"));
                        pb.inc(1);
                    }
                }
            }
        }
        Err(e) => {
            errors
                .lock()
                .unwrap()
                .push(format!("Error reading directory entry: {e}"));
        }
    });

    pb.finish_with_message("Processing complete!");

    // Report errors at the end
    let errors = errors.into_inner().unwrap();
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

fn process_single_file(
    input_file: String,
    output_file: String,
    output_depth: SASALevel,
    n_points: usize,
) -> Result<(), CLIError> {
    println!("Processing single file...");

    let format = OutputFormat::from_file_extension(&output_file);

    process(
        input_file,
        output_file,
        output_depth,
        &format,
        n_points,
        true,
    )?;
    println!("Finished!");
    Ok(())
}

fn main() -> Result<(), CLIError> {
    let args = Args::parse();

    if std::path::Path::new(&args.input).is_dir() {
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
            &format,
        )
    } else {
        process_single_file(args.input, args.output, args.output_depth, args.n_points)
    }
}
