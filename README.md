# ampSeqPrimerDesign

A Python script for designing amplicon sequencing primer sets at user-defined genomic target sites. Primers are designed using the Primer3 engine and automatically adapted with SP5/SP7 Illumina-compatible handle sequences. Output is written to two CSV files — one containing the bare gene-specific primers (GSPs) and one containing the fully adapted sequences ready for ordering.

---

## Overview

The script accepts a standard BED file of target coordinates, expands each interval into a design space, fetches the corresponding reference genome sequence using pybedtools, and submits each sequence to Primer3 for primer design. A deduplication step ensures that all reported primer sets for a given locus are sequence-unique, giving you real choices rather than redundant candidates. Final output includes both the raw GSP sequences and handle-appended sequences in separate CSVs.

---

## Requirements

- Python 3.10
- Conda (recommended for environment management)
- `bedtools` accessible on `PATH`
- A BWA-indexed reference genome FASTA (e.g. hg38)

### Conda Environment

An `environment.yaml` is written automatically to your output directory on each run. You can also create it manually before first use:

```bash
conda env create -f environment.yaml
conda activate ampseq_primer_design
```

The environment includes: `primer3-py`, `pandas`, `pybedtools`, `bedtools`, `pyyaml`

---

## Input Format

The script expects a **tab-delimited BED file with no header** and the following six columns:

| Column | Description |
|--------|-------------|
| `chrom` | Chromosome name (e.g. `chr2`) |
| `start` | 0-based start coordinate of the target site |
| `stop` | End coordinate of the target site |
| `name` | Locus name/identifier (used in output labeling) |
| `score` | BED score field (unused, but required for format compliance) |
| `strand` | Strand (`+` or `-`, unused by the script) |

Example:
```
chr2    73160351    73160371    CCR5_Site1    0    +
chr17   41244000    41244020    BRCA1_Site1   0    +
```

> The target coordinates should represent the center of your desired amplicon (e.g. the cut site). The script automatically expands each interval by the design space buffer to define the primer search region.

---

## Usage

```bash
python ampSeqPrimerDesign.py \
  -f targets.bed \
  -r /path/to/hg38.fa \
  -o Outputs/ \
  -s 200 \
  -b 50 \
  -n 3
```

### All Parameters

| Flag | Long form | Default | Description |
|------|-----------|---------|-------------|
| `-f` | `--fileName` | *(required)* | Input BED file of target sites |
| `-r` | `--referencePath` | *(required)* | Path to reference genome FASTA (must be accessible by pybedtools) |
| `-o` | `--outDir` | `Outputs` | Directory for all output files and working intermediates |
| `-s` | `--ampliconSize` | `200` | Target amplicon size in base pairs |
| `-b` | `--buffer` | `50` | Acceptable size buffer around target amplicon; Primer3 will accept products in the range `[ampliconSize - buffer, ampliconSize + buffer]` |
| `-d` | `--designSpace` | `1.5` | Ratio of design space to amplicon size. Must be `>1.0`. A value of `1.5` on a 200 bp amplicon gives a 300 bp search window (150 bp on each side of the target). |
| `-t` | `--meltingTemperature` | `60` | Optimal primer Tm in °C. Min/max Tm are set to ±3°C of this value. |
| `-n` | `--numberSets` | `3` | Number of unique primer sets to return per locus |
| `-p` | `--probe01` | `0` | Set to `1` to also design an internal probe/oligo for each locus |
| `--FWDHandle` | | `ACACTCTTTCCCTACACGACGCTCTTCCGATCT` | SP5 handle sequence prepended to forward primers |
| `--REVHandle` | | `GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT` | SP7 handle sequence prepended to reverse primers |

---

## How It Works

### 1. Design Space Expansion
Each target interval is expanded symmetrically by `(ampliconSize × designSpace) / 2` base pairs on each side. These expanded coordinates are written to a temporary BED file and passed to `pybedtools` to fetch the corresponding FASTA sequence from the reference genome.

### 2. Primer Design (Primer3)
Each FASTA sequence is submitted to Primer3 via `primer3-py`. To support deduplication, the script requests a candidate pool of `min(numberSets × 5, 25)` primer pairs internally, even if fewer unique sets are ultimately reported. Key Primer3 parameters:

| Parameter | Value |
|-----------|-------|
| Optimal primer length | 20 bp |
| Primer length range | 17–25 bp |
| Optimal Tm | user-defined |
| Tm range | ±3°C of optimal |
| Probe Tm | optimal + 5°C (±3°C) |
| GC content | 20–80% |
| Max poly-X run | 5 |
| Salt (monovalent) | 50 mM |
| DNA concentration | 50 nM |

### 3. Deduplication
Primer3 ranked candidates are walked in order of penalty score. A candidate pair is accepted only if neither its forward nor reverse sequence has been seen before for that locus. This continues until `numberSets` unique pairs are collected or the candidate pool is exhausted. A warning is printed if fewer unique sets than requested are available.

### 4. Handle Appending
Accepted GSP sequences are copied to a second dataframe and prepended with the SP5 (FWD) and SP7 (REV) handle sequences to produce the full adapted primers ready for synthesis ordering.

---

## Output Files

Each run writes all output into an automatically created subdirectory of `--outDir` structured as `<outDir>/YYYYMMDD/RunX`, where `YYYYMMDD` is today's date and `X` increments based on how many runs have already been completed that day. For example, the third run on April 26 2026 against `Outputs/` would produce:

```
Outputs/
└── 20260426/
    ├── Run1/
    ├── Run2/
    └── Run3/
```

The resolved run directory is printed to the console at the start of each run.

| File | Description |
|------|-------------|
| `DesignedPrimers_GSP_<timestamp>.csv` | Gene-specific primer sequences only, with Tm and product metrics |
| `DesignedPrimers_Adapted_<timestamp>.csv` | Full primers with SP5/SP7 handles prepended |
| `DesignSpaceBED_<timestamp>.txt` | Expanded BED intervals used for sequence retrieval (intermediate) |
| `DesignSequences_<timestamp>.txt` | FASTA sequences fetched from reference (intermediate) |
| `environment.yaml` | Conda environment file for reproducing the script's dependencies |

### CSV Column Descriptions

| Column | Description |
|--------|-------------|
| `Site` | Locus name from the input BED `name` field |
| `Set Number` | Primer set index (Set 1, Set 2, ...) ordered by Primer3 penalty score |
| `Left Primer` | Forward (left) primer sequence |
| `Left Primer Tm (C)` | Predicted Tm of the forward primer, rounded to the nearest 0.1°C |
| `Right Primer` | Reverse (right) primer sequence |
| `Right Primer Tm (C)` | Predicted Tm of the reverse primer, rounded to the nearest 0.1°C |
| `Probe` | Internal probe sequence (if `-p 1` was set, otherwise `No Probe Designed`) |
| `Probe Tm (C)` | Predicted probe Tm, rounded to the nearest 0.1°C |
| `Product Size` | Expected amplicon size in base pairs |
| `Product Penalty Score` | Primer3 pair penalty score, rounded to the nearest 0.001 (lower is better) |

> In the **Adapted** CSV, the `Left Primer` and `Right Primer` columns contain the full handle + GSP sequence. The GSP portion begins after the handle (33 nt for SP5, 34 nt for SP7 using the default handles).

---

## Notes and Caveats

- **bedtools must be on PATH.** The script calls `shutil.which('bedtools')` at startup and exits with a clear error if it is not found. The recommended approach is to run the script from within the activated conda environment.
- **Reference genome indexing.** pybedtools uses `bedtools getfasta` internally, which requires the FASTA to be accessible but does not require a pre-built BWA index. A `.fai` samtools index on the FASTA is recommended for performance.
- **Design space ratio.** If `--designSpace` is too close to `1.0`, the primer search window may be too narrow for Primer3 to find valid pairs at the requested amplicon size. A value between `1.3` and `2.0` is generally appropriate.
- **Duplicate warnings.** If a locus receives a warning that fewer unique sets were found than requested, consider increasing `--designSpace` to give Primer3 a larger search window, or relaxing the Tm or size constraints.
- **Probe design.** When `-p 1` is set, Primer3 will attempt to design an internal oligo at an optimum Tm of `meltingTemperature + 5°C`. This is suitable for hydrolysis probe assays. Probe sequences are included in both output CSVs.
