use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_bed_nuc::{NucOptions, nuc};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

/// bedtools nuc uses multi-char single-dash flags; clap requires long flags, so:
/// `-fi`→`--fasta`, `-bed`→`--bed`, `-s`→`--strand`, `-seq`→`--seq`,
/// `-C`→`--case-insensitive`, `-fullHeader`→`--full-header`, `-pattern`→`--pattern`.
#[allow(clippy::struct_excessive_bools)] // five orthogonal boolean flags; no state machine applies
#[derive(Parser, Debug)]
#[command(
    name = "rsomics-bed-nuc",
    version,
    about = "Per-interval nucleotide composition from FASTA + BED",
    long_about = None,
    disable_help_flag = true
)]
pub struct Cli {
    /// Indexed FASTA file (.fai must exist alongside it; run `samtools faidx`).
    #[arg(long = "fasta", value_name = "FILE")]
    pub fasta: PathBuf,

    /// BED file of intervals to profile.
    #[arg(long = "bed", value_name = "FILE")]
    pub bed: PathBuf,

    /// Reverse-complement minus-strand intervals before counting (requires BED col 6).
    #[arg(long = "strand", short = 's')]
    pub strand: bool,

    /// Append the extracted reference sequence to each output line.
    #[arg(long = "seq")]
    pub seq: bool,

    /// Use the full FASTA header (not just first word) for sequence lookup.
    #[arg(long = "full-header")]
    pub full_header: bool,

    /// Count occurrences of this exact pattern in each interval.
    #[arg(long = "pattern", value_name = "PATTERN")]
    pub pattern: Option<String>,

    /// Match --pattern case-insensitively.
    #[arg(long = "case-insensitive", short = 'C')]
    pub case_insensitive: bool,

    /// Output file (default: stdout).
    #[arg(short = 'o', long, default_value = "-")]
    pub output: String,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }

    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let opts = NucOptions {
            strand_aware: self.strand,
            print_seq: self.seq,
            full_header: self.full_header,
            pattern: self.pattern.as_deref(),
            case_insensitive: self.case_insensitive,
        };

        let mut out: Box<dyn std::io::Write> = if self.output == "-" {
            Box::new(std::io::stdout().lock())
        } else {
            Box::new(std::fs::File::create(&self.output).map_err(rsomics_common::RsomicsError::Io)?)
        };

        nuc(&self.fasta, &self.bed, &opts, &mut out)
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
    tagline: "Per-interval nucleotide composition from indexed FASTA + BED (bedtools nuc equivalent).",
    origin: Some(Origin {
        upstream: "bedtools nuc",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &[
        "--fasta <ref.fa> --bed <regions.bed> [-s] [--seq] [--pattern PAT] [-o out.tsv]",
    ],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: None,
                long: "fasta",
                aliases: &["-fi"],
                value: Some("<path>"),
                type_hint: Some("PathBuf"),
                required: true,
                default: None,
                description: "Indexed FASTA file (.fai must exist alongside it).",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "bed",
                aliases: &["-bed"],
                value: Some("<path>"),
                type_hint: Some("PathBuf"),
                required: true,
                default: None,
                description: "BED intervals to profile.",
                why_default: None,
            },
            FlagSpec {
                short: Some('s'),
                long: "strand",
                aliases: &["-s"],
                value: None,
                type_hint: None,
                required: false,
                default: None,
                description: "Reverse-complement minus-strand intervals before counting (requires BED col 6).",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "seq",
                aliases: &["-seq"],
                value: None,
                type_hint: None,
                required: false,
                default: None,
                description: "Append the extracted reference sequence to each output row.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "full-header",
                aliases: &["-fullHeader"],
                value: None,
                type_hint: None,
                required: false,
                default: None,
                description: "Use full FASTA header (not just first word) for sequence lookup.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "pattern",
                aliases: &["-pattern"],
                value: Some("PAT"),
                type_hint: Some("String"),
                required: false,
                default: None,
                description: "Count occurrences of this exact pattern per interval.",
                why_default: None,
            },
            FlagSpec {
                short: Some('C'),
                long: "case-insensitive",
                aliases: &["-C"],
                value: None,
                type_hint: None,
                required: false,
                default: None,
                description: "Match --pattern case-insensitively.",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "GC content of ATAC-seq peaks",
            command: "rsomics-bed-nuc --fasta hg38.fa --bed peaks.bed",
        },
        Example {
            description: "Strand-aware composition with extracted sequence",
            command: "rsomics-bed-nuc --fasta hg38.fa --bed motifs.bed --strand --seq",
        },
        Example {
            description: "Count CpG dinucleotides in promoter windows",
            command: "rsomics-bed-nuc --fasta hg38.fa --bed promoters.bed --pattern CG",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
