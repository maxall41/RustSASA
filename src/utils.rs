use snafu::Snafu;

const LANES: usize = 16;
pub fn simd_sum(values: &[f32]) -> f32 {
    let chunks = values.chunks_exact(LANES);
    let remainder = chunks.remainder();

    let sum = chunks.fold([0.0f32; LANES], |mut acc, chunk| {
        let chunk: [f32; LANES] = chunk.try_into().unwrap();
        for i in 0..LANES {
            acc[i] += chunk[i];
        }
        acc
    });

    let remainder: f32 = remainder.iter().copied().sum();

    let mut reduced = 0.0f32;
    for i in 0..LANES {
        reduced += sum[i];
    }
    reduced + remainder
}

pub(crate) fn serialize_chain_id(s: &str) -> isize {
    let mut result = 0;
    for c in s.chars() {
        if c.is_ascii_alphabetic() {
            let position = c.to_ascii_uppercase() as isize - 64;
            result = result * 10 + position;
        }
    }
    result
}