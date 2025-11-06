pub mod consts;
pub mod io;
use std::sync::LazyLock;

use pulp::Arch;

static ARCH: LazyLock<Arch> = LazyLock::new(Arch::new);

pub fn simd_sum(values: &[f32]) -> f32 {
    let mut total = 0f32;
    ARCH.dispatch(|| {
        for x in values {
            total += x;
        }
    });
    total
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
