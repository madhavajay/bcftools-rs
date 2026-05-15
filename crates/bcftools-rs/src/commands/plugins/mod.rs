//! In-process plugin record-processing implementations.
//!
//! Each module mirrors one upstream `bcftools/plugins/<name>.c` and is invoked
//! through `bcftools plugin <name>` / `bcftools +<name>`. The static
//! registry/listing surface lives in [`super::plugin`]; this module holds the
//! actual per-plugin algorithms as they are ported.

pub mod counts;
pub mod fill_an_ac;
pub mod missing2ref;
