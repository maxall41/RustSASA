use rust_sasa::SASAResult;
use std::fs;
use std::path::Path;

pub fn read_json_result(path: &Path) -> SASAResult {
    let content = fs::read_to_string(path).expect("Failed to read output file");
    serde_json::from_str(&content).expect("Failed to deserialize JSON content")
}
