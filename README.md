# rsomics-bed-nuc

Profile nucleotide content of BED intervals in a FASTA file — a fast Rust
reimplementation of `bedtools nuc`.

## Usage

```
rsomics-bed-nuc --fi genome.fa --bed regions.bed
```

Options:
- `--fi <fasta>` — indexed FASTA file (requires `.fai` sidecar)
- `--bed <bed>` — input BED file (default: stdin)
- `-s / --strand` — profile according to strand (reverse-complement `-` strand intervals)
- `--seq` — print the extracted sequence as an extra column
- `--pattern <str>` — count occurrences of this pattern in each interval
- `-C` — case-insensitive pattern matching

## Install

```
cargo install rsomics-bed-nuc
```

## Origin

This crate is an independent Rust reimplementation of `bedtools nuc` based on:
- Quinlan & Hall (2010). BEDTools: a flexible suite of utilities for comparing
  genomic features. Bioinformatics 26(6): 841–842. DOI: 10.1093/bioinformatics/btq033
- The FASTA/BED format specifications and black-box behavior testing against
  `bedtools nuc 2.31.1`

No source code from the BEDTools upstream was used as reference during implementation.
Test fixtures are independently generated.

License: MIT OR Apache-2.0.
Upstream credit: BEDTools <https://github.com/arq5x/bedtools2> (GPL-2.0).
