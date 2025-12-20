// Copyright (c) 2024 Maxwell Campbell. Licensed under the MIT License.
use rust_sasa::SASAResult;
use serde::Deserialize;
use std::fs;
use std::path::Path;

#[allow(dead_code)] // Not actually dead code. Rust does not recognize its usage in tests for some reason.
pub fn read_json_result(path: &Path) -> SASAResult {
    let content = fs::read_to_string(path).expect("Failed to read output file");
    serde_json::from_str(&content).expect("Failed to deserialize JSON content")
}

#[allow(dead_code)]
pub fn read_xml_result(path: &Path) -> SASAResult {
    let content = fs::read_to_string(path).expect("Failed to read output file");
    quick_xml::de::from_str(&content).expect("Failed to deserialize XML content")
}

// FreeSASA JSON structures
#[derive(Debug, Deserialize)]
#[allow(dead_code)] // Not actually dead code. Rust does not recognize its usage in tests for some reason.
pub struct FreeSASAOutput {
    pub results: Vec<FreeSASAResult>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)] // Not actually dead code. Rust does not recognize its usage in tests for some reason.
pub struct FreeSASAResult {
    pub structure: Vec<FreeSASAStructure>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)] // Not actually dead code. Rust does not recognize its usage in tests for some reason.
pub struct FreeSASAStructure {
    pub chains: Vec<FreeSASAChain>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)] // Not actually dead code. Rust does not recognize its usage in tests for some reason.
pub struct FreeSASAChain {
    pub label: String,
    pub area: FreeSASAArea,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)] // Not actually dead code. Rust does not recognize its usage in tests for some reason.
pub struct FreeSASAArea {
    pub total: f64,
}
