//! Custom benchmark smoke harness for bcftools-rs command work.
//!
//! Run with:
//!
//! ```text
//! cargo bench -p bcftools-rs --bench command_harness
//! ```
//!
//! The harness keeps named slots for the parity-critical areas tracked in
//! `TODO.md`: view, merge, norm, annotate, and the filter engine. Slots whose
//! implementation is not present yet report as blocked instead of pretending to
//! have meaningful performance data.

use std::fs::File;
use std::hint::black_box;
use std::io::BufReader;
use std::path::Path;
use std::time::{Duration, Instant};

use tempfile::TempDir;

const ITERATIONS: usize = 20;

const VCF: &str = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t10\t.\tA\tG\t50\tPASS\t.\n\
1\t20\t.\tC\tT\t60\tPASS\t.\n\
1\t30\t.\tG\tA\t70\tPASS\t.\n";

fn main() {
    let temp = TempDir::new().expect("tempdir");
    let input = temp.path().join("input.vcf");
    std::fs::write(&input, VCF).expect("write fixture");

    println!("bcftools-rs custom command bench");
    bench_view_read_path(&input);
    blocked(
        "merge",
        "command port depends on synced reader and bcf_translate coverage",
    );
    blocked(
        "norm",
        "command port depends on abuf, vcfbuf, and reference-check logic",
    );
    blocked(
        "annotate",
        "command port depends on synced reader, regidx, and bcf_update_* coverage",
    );
    blocked(
        "filter-engine",
        "expression compiler/evaluator is not ported yet",
    );
}

fn bench_view_read_path(input: &Path) {
    let elapsed = bench("view_read_path", || {
        let format = htslib_rs::format::detect_path(input).expect("detect VCF");
        black_box(format);

        let file = File::open(input).expect("open VCF");
        let mut reader = htslib_rs::vcf::io::Reader::new(BufReader::new(file));
        let header = reader.read_header().expect("read VCF header");
        let count = reader.records().count();
        black_box((header, count));
    });

    report("view", elapsed);
}

fn bench<F>(name: &str, mut f: F) -> Duration
where
    F: FnMut(),
{
    let start = Instant::now();
    for _ in 0..ITERATIONS {
        f();
    }
    let elapsed = start.elapsed();
    black_box(name);
    elapsed / ITERATIONS as u32
}

fn report(name: &str, mean: Duration) {
    println!("{name:24} mean={mean:?} iterations={ITERATIONS}");
}

fn blocked(name: &str, reason: &str) {
    println!("{name:24} blocked: {reason}");
}
