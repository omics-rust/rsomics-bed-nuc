use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_nuc(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-bed-nuc");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let fasta = manifest.join("tests/golden/ref.fa");
    let bed = manifest.join("tests/golden/regions.bed");

    c.bench_function("rsomics-bed-nuc golden", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .args([
                    "--fi",
                    fasta.to_str().unwrap(),
                    "--bed",
                    bed.to_str().unwrap(),
                ])
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_nuc);
criterion_main!(benches);
