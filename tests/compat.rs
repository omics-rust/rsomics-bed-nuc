//! Compatibility tests: byte-identical output vs bedtools nuc v2.31.1.
//!
//! Golden fixtures generated with:
//!   bedtools nuc -fi tests/golden/ref.fa -bed tests/golden/regions.bed
//!   bedtools nuc -fi tests/golden/ref.fa -bed tests/golden/strand.bed -s
//!   bedtools nuc -fi tests/golden/ref.fa -bed tests/golden/regions.bed -seq
//!   bedtools nuc -fi tests/golden/ref.fa -bed tests/golden/regions.bed -pattern ACGT
//!   bedtools nuc -fi tests/golden/ref.fa -bed tests/golden/regions.bed -pattern acgt -C

use std::path::Path;

use rsomics_bed_nuc::{NucOptions, nuc};

fn golden(name: &str) -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn run(fasta: &str, bed: &str, opts: &NucOptions<'_>) -> String {
    let mut buf = Vec::new();
    nuc(&golden(fasta), &golden(bed), opts, &mut buf).expect("nuc failed");
    String::from_utf8(buf).expect("valid utf8")
}

fn expected(name: &str) -> String {
    std::fs::read_to_string(golden(name)).expect("golden file missing")
}

fn basic_opts() -> NucOptions<'static> {
    NucOptions {
        strand_aware: false,
        print_seq: false,
        full_header: false,
        pattern: None,
        case_insensitive: false,
    }
}

#[test]
fn basic_composition() {
    let opts = basic_opts();
    assert_eq!(
        run("ref.fa", "regions.bed", &opts),
        expected("basic.tsv"),
        "basic composition mismatch"
    );
}

#[test]
fn strand_aware() {
    let opts = NucOptions {
        strand_aware: true,
        ..basic_opts()
    };
    assert_eq!(
        run("ref.fa", "strand.bed", &opts),
        expected("strand_aware.tsv"),
        "strand-aware mismatch"
    );
}

#[test]
fn with_sequence() {
    let opts = NucOptions {
        print_seq: true,
        ..basic_opts()
    };
    assert_eq!(
        run("ref.fa", "regions.bed", &opts),
        expected("with_seq.tsv"),
        "-seq output mismatch"
    );
}

#[test]
fn pattern_case_sensitive() {
    let opts = NucOptions {
        pattern: Some("ACGT"),
        ..basic_opts()
    };
    assert_eq!(
        run("ref.fa", "regions.bed", &opts),
        expected("with_pattern.tsv"),
        "-pattern (case-sensitive) mismatch"
    );
}

#[test]
fn pattern_case_insensitive() {
    let opts = NucOptions {
        pattern: Some("acgt"),
        case_insensitive: true,
        ..basic_opts()
    };
    assert_eq!(
        run("ref.fa", "regions.bed", &opts),
        expected("with_pattern_ci.tsv"),
        "-pattern -C (case-insensitive) mismatch"
    );
}
