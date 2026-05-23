//! Per-interval nucleotide composition from FASTA + BED.
//!
//! Algorithm: load the indexed FASTA entirely into memory (one `Vec<u8>` per
//! sequence, newlines stripped), then scan each BED interval in O(L) where L
//! is the interval length. Total complexity is O(G + N·L̄) where G = genome
//! size, N = interval count, L̄ = mean interval length.  No per-element
//! rescan; no O(N²) pathology.
//!
//! Output mirrors `bedtools nuc`:
//!   <all BED columns...>\t%AT\t%GC\t#A\t#C\t#G\t#T\t#N\t#other\tlen
//! With optional extra columns for --seq and --pattern.

#![allow(clippy::cast_precision_loss)] // u64→f64 for pct_at/pct_gc intentional

use std::collections::HashMap;
use std::fmt::Write as FmtWrite;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

/// Options for `nuc`.
#[allow(clippy::struct_excessive_bools)] // five orthogonal flags, no state machine applies
pub struct NucOptions<'a> {
    /// Profile sequence according to strand (reverse-complement minus strand).
    pub strand_aware: bool,
    /// Append the extracted sequence to each output line.
    pub print_seq: bool,
    /// Use the full FASTA header (not just the first word).
    pub full_header: bool,
    /// Count occurrences of this exact pattern in the interval sequence.
    pub pattern: Option<&'a str>,
    /// Match pattern case-insensitively (only effective when pattern is set).
    pub case_insensitive: bool,
}

/// Per-interval nucleotide statistics.
pub struct NucStats {
    pub pct_at: f64,
    pub pct_gc: f64,
    pub num_a: u64,
    pub num_c: u64,
    pub num_g: u64,
    pub num_t: u64,
    pub num_n: u64,
    pub num_other: u64,
    pub seq_len: u64,
    /// Pattern match count (only populated when a pattern is specified).
    pub pattern_count: Option<u64>,
}

/// FAI entry: just enough to seek into the FASTA file.
struct FaiEntry {
    name: String,
    length: u64,
    offset: u64,
    line_bases: u64,
    line_width: u64,
}

/// Load a `.fai` index file, returning entries in file order.
fn read_fai(fai_path: &Path) -> Result<Vec<FaiEntry>> {
    let file = File::open(fai_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", fai_path.display())))?;
    let reader = BufReader::new(file);
    let mut entries = Vec::new();
    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 5 {
            return Err(RsomicsError::InvalidInput(format!(
                "malformed .fai line (need 5 fields): {line}"
            )));
        }
        let parse = |s: &str, field: &str| -> Result<u64> {
            s.parse()
                .map_err(|e| RsomicsError::InvalidInput(format!(".fai {field} parse error: {e}")))
        };
        entries.push(FaiEntry {
            name: cols[0].to_string(),
            length: parse(cols[1], "length")?,
            offset: parse(cols[2], "offset")?,
            line_bases: parse(cols[3], "line_bases")?,
            line_width: parse(cols[4], "line_width")?,
        });
    }
    Ok(entries)
}

/// Load the entire FASTA into a `HashMap<name -> Vec<u8>>` using the FAI for
/// fast seeking. Newlines are stripped; bases stored as-is (mixed case).
fn load_fasta(fasta_path: &Path, fai: &[FaiEntry]) -> Result<HashMap<String, Vec<u8>>> {
    let mut file = File::open(fasta_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", fasta_path.display())))?;
    let mut map: HashMap<String, Vec<u8>> = HashMap::with_capacity(fai.len());

    for entry in fai {
        #[allow(clippy::cast_possible_truncation)]
        let mut seq: Vec<u8> = Vec::with_capacity(entry.length as usize);
        file.seek(SeekFrom::Start(entry.offset))
            .map_err(RsomicsError::Io)?;

        // Read enough raw bytes to span all sequence lines for this entry.
        // line_width includes the newline; total raw bytes ≤ ceil(length/line_bases)*line_width.
        let lines = entry.length.div_ceil(entry.line_bases);
        let raw_bytes = lines * entry.line_width;
        #[allow(clippy::cast_possible_truncation)]
        let mut buf = vec![0u8; raw_bytes as usize];
        let n = file.read(&mut buf).map_err(RsomicsError::Io)?;
        for &b in &buf[..n] {
            if b != b'\n' && b != b'\r' {
                seq.push(b);
                if seq.len() as u64 == entry.length {
                    break;
                }
            }
        }

        map.insert(entry.name.clone(), seq);
    }
    Ok(map)
}

/// Reverse-complement a sequence slice in place.
fn revcomp(seq: &mut [u8]) {
    seq.reverse();
    for b in seq.iter_mut() {
        *b = match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            other => other, // N and ambiguity codes stay as uppercase
        };
    }
}

/// Count overlapping occurrences of `pattern` in `seq` (Knuth-Morris-Pratt).
fn count_pattern(seq: &[u8], pattern: &str, case_insensitive: bool) -> u64 {
    let pat: Vec<u8> = if case_insensitive {
        pattern.bytes().map(|b| b.to_ascii_uppercase()).collect()
    } else {
        pattern.bytes().collect()
    };
    if pat.is_empty() || seq.len() < pat.len() {
        return 0;
    }
    // Build KMP failure table.
    let m = pat.len();
    let mut fail = vec![0usize; m];
    let mut k = 0usize;
    for i in 1..m {
        while k > 0 && pat[k] != pat[i] {
            k = fail[k - 1];
        }
        if pat[k] == pat[i] {
            k += 1;
        }
        fail[i] = k;
    }
    // Search, allowing overlapping matches.
    let mut count = 0u64;
    let mut q = 0usize;
    for &b in seq {
        let b = if case_insensitive {
            b.to_ascii_uppercase()
        } else {
            b
        };
        while q > 0 && pat[q] != b {
            q = fail[q - 1];
        }
        if pat[q] == b {
            q += 1;
        }
        if q == m {
            count += 1;
            q = fail[q - 1]; // allow overlapping matches
        }
    }
    count
}

/// Compute nucleotide stats for a single interval sequence slice.
#[must_use]
pub fn compute_stats(bases: &[u8], opts: &NucOptions<'_>) -> NucStats {
    let mut num_a = 0u64;
    let mut num_c = 0u64;
    let mut num_g = 0u64;
    let mut num_t = 0u64;
    let mut num_n = 0u64;
    let mut num_other = 0u64;
    for &b in bases {
        match b.to_ascii_uppercase() {
            b'A' => num_a += 1,
            b'C' => num_c += 1,
            b'G' => num_g += 1,
            b'T' => num_t += 1,
            b'N' => num_n += 1,
            _ => num_other += 1,
        }
    }
    let seq_len = bases.len() as u64;
    // bedtools nuc computes: printf("%f\t%f\t", (float)(a+t)/seqLength, (float)(c+g)/seqLength)
    // where seqLength is int64_t. In C, float / int64_t promotes the int64_t to float, so
    // the entire division is done in f32 hardware. printf's %f then promotes the f32 result
    // to f64 for printing. We replicate this: cast numerator to f32, cast denominator to f32,
    // divide in f32, then widen to f64 for the {:.6} format — matching byte-for-byte.
    let (pct_at, pct_gc) = if seq_len == 0 {
        (0.0f64, 0.0f64)
    } else {
        let len_f32 = seq_len as f32;
        let at = f64::from((num_a + num_t) as f32 / len_f32);
        let gc = f64::from((num_c + num_g) as f32 / len_f32);
        (at, gc)
    };
    let pattern_count = opts
        .pattern
        .map(|p| count_pattern(bases, p, opts.case_insensitive));
    NucStats {
        pct_at,
        pct_gc,
        num_a,
        num_c,
        num_g,
        num_t,
        num_n,
        num_other,
        seq_len,
        pattern_count,
    }
}

/// Write the header line to `out`.
///
/// bedtools nuc prints `#1_usercol\t2_usercol\t...\tN_pct_at\t...`.
/// Column numbering is determined by the BED column count of the first data line.
fn write_header(num_bed_cols: usize, opts: &NucOptions<'_>, out: &mut impl Write) -> Result<()> {
    let col_offset = num_bed_cols + 1;
    let mut header = String::new();
    for i in 1..=num_bed_cols {
        if i > 1 {
            header.push('\t');
        }
        write!(header, "{i}_usercol").unwrap();
    }
    let at = col_offset;
    let gc = at + 1;
    let na = gc + 1;
    let nc = na + 1;
    let ng = nc + 1;
    let nt = ng + 1;
    let nn = nt + 1;
    let no = nn + 1;
    let nl = no + 1;
    write!(
        header,
        "\t{at}_pct_at\t{gc}_pct_gc\t{na}_num_A\t{nc}_num_C\t{ng}_num_G\t{nt}_num_T\t{nn}_num_N\t{no}_num_oth\t{nl}_seq_len"
    ).unwrap();
    if opts.print_seq {
        let ns = nl + 1;
        write!(header, "\t{ns}_seq").unwrap();
    }
    if opts.pattern.is_some() {
        let np = nl + if opts.print_seq { 2 } else { 1 };
        write!(header, "\t{np}_user_patt_count").unwrap();
    }
    writeln!(out, "#{header}").map_err(RsomicsError::Io)
}

/// Run `bedtools nuc` equivalent: for each BED interval, extract the FASTA
/// subsequence and emit nucleotide composition columns.
pub fn nuc(
    fasta_path: &Path,
    bed_path: &Path,
    opts: &NucOptions<'_>,
    out: &mut dyn Write,
) -> Result<()> {
    // Require a .fai index alongside the FASTA (samtools faidx produces it).
    let fai_path = {
        let mut p = fasta_path.as_os_str().to_os_string();
        p.push(".fai");
        std::path::PathBuf::from(p)
    };
    let fai = read_fai(&fai_path)?;
    let genome = load_fasta(fasta_path, &fai)?;

    // Length map for boundary checking (mirrors bedtools nuc stderr skip logic).
    let lengths: HashMap<&str, u64> = fai.iter().map(|e| (e.name.as_str(), e.length)).collect();

    let bed_file = File::open(bed_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", bed_path.display())))?;
    let reader = BufReader::new(bed_file);
    let mut out = BufWriter::with_capacity(256 * 1024, out);

    let mut header_written = false;

    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        let line = line.trim_end_matches('\r');
        if line.is_empty()
            || line.starts_with('#')
            || line.starts_with("track")
            || line.starts_with("browser")
        {
            continue;
        }

        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 3 {
            continue;
        }

        // Emit header using the column count of the first data line.
        if !header_written {
            write_header(cols.len(), opts, &mut out)?;
            header_written = true;
        }

        let chrom = cols[0];
        let start: u64 = cols[1]
            .parse()
            .map_err(|e| RsomicsError::InvalidInput(format!("BED start parse: {e}")))?;
        let end: u64 = cols[2]
            .parse()
            .map_err(|e| RsomicsError::InvalidInput(format!("BED end parse: {e}")))?;

        // bedtools nuc skips (with a stderr message) intervals that exceed the chromosome length.
        let chrom_len = lengths.get(chrom).copied().unwrap_or(0);
        if end > chrom_len {
            eprintln!(
                "Feature ({chrom}:{start}-{end}) beyond the length of {chrom} size ({chrom_len} bp).  Skipping."
            );
            continue;
        }

        let Some(seq_full) = genome.get(chrom) else {
            eprintln!("WARNING: sequence {chrom} not found in FASTA index. Skipping.");
            continue;
        };

        #[allow(clippy::cast_possible_truncation)]
        let (s, e) = (start as usize, end as usize);
        let mut bases: Vec<u8> = seq_full[s..e].to_vec();

        // Reverse-complement minus-strand intervals when -s is active.
        if opts.strand_aware {
            let strand = if cols.len() >= 6 { cols[5] } else { "." };
            if strand == "-" {
                revcomp(&mut bases);
            }
        }

        let stats = compute_stats(&bases, opts);

        // Emit all original BED columns followed by the stat columns.
        for (i, col) in cols.iter().enumerate() {
            if i > 0 {
                out.write_all(b"\t").map_err(RsomicsError::Io)?;
            }
            out.write_all(col.as_bytes()).map_err(RsomicsError::Io)?;
        }
        write!(
            out,
            "\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            stats.pct_at,
            stats.pct_gc,
            stats.num_a,
            stats.num_c,
            stats.num_g,
            stats.num_t,
            stats.num_n,
            stats.num_other,
            stats.seq_len,
        )
        .map_err(RsomicsError::Io)?;
        if opts.print_seq {
            // Output the sequence as it was processed (post revcomp if -s active).
            out.write_all(b"\t").map_err(RsomicsError::Io)?;
            out.write_all(&bases).map_err(RsomicsError::Io)?;
        }
        if let Some(pc) = stats.pattern_count {
            write!(out, "\t{pc}").map_err(RsomicsError::Io)?;
        }
        out.write_all(b"\n").map_err(RsomicsError::Io)?;
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(())
}
