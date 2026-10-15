# CARDOAP

CARDOAP integrates the CARD protein homolog database with the ARGs-OAP framework for short-read ARG profiling and cell-normalized quantification.

The current database contains **6,059 CARD protein homolog reference sequences**. CARDOAP retains the original two-stage ARGs-OAP workflow and KO30-based cell normalization, while replacing the default SARG database with CARD. CARD references associated with multiple drug classes were individually checked and validated to ensure consistent class assignment and to avoid duplicate counting during quantification. The database also includes **128 broad-spectrum multidrug efflux resistance genes**, which were specifically reviewed and classified as multidrug resistance determinants.

## Installation

First install ARGs-OAP following the official instructions:

https://github.com/xinehc/args_oap

Then clone CARDOAP:

```bash
git clone https://github.com/uby76/CARDOAP.git
cd CARDOAP
```

## Build the CARDOAP database

The CARD protein FASTA is provided in:

```text
database/CARD_protein_homolog_full_6059.fasta
```

Before the first run, build the DIAMOND and NCBI BLAST databases:

```bash
args_oap make_db \
  -i database/CARD_protein_homolog_full_6059.fasta
```

This step only needs to be performed once.

## Quick test

A small paired-end test dataset is provided in the `test/` directory.

Run the complete test workflow:

```bash
bash scripts/run_test.sh
```

The script automatically performs database indexing, Stage 1 screening, Stage 2 ARG identification, and abundance normalization.

## Run CARDOAP on your own reads

Prepare paired-end reads using names such as:

```text
sample1_1.fa
sample1_2.fa
sample2_1.fa
sample2_2.fa
```

Place the input files in an input directory, for example:

```text
input/
├── sample1_1.fa
├── sample1_2.fa
├── sample2_1.fa
└── sample2_2.fa
```

### Stage 1

Run DIAMOND BLASTX screening against the CARD database:

```bash
args_oap stage_one \
  -i input \
  -o output \
  -f fa \
  -t 8 \
  --database database/CARD_protein_homolog_full_6059.fasta
```

Stage 1 identifies candidate ARG-like reads and estimates sample-level normalization factors including `n16S` and `nCell`.

### Stage 2

Run stringent NCBI BLASTX identification and quantification:

```bash
args_oap stage_two \
  -i output \
  -t 8 \
  --database database/CARD_protein_homolog_full_6059.fasta \
  --structure1 database/CARD_structure_full_6059.txt
```

## Main output

For ARG copies per cell, use:

```text
output/normalized_cell.reference.txt
output/normalized_cell.aro.txt
output/normalized_cell.arg.txt
output/normalized_cell.family.txt
output/normalized_cell.class.txt
output/normalized_cell.mechanism.txt
```

For example:

```text
normalized_cell.arg.txt
```

contains CARD ARG abundances expressed as cell-normalized copy signals using the ARGs-OAP KO30 normalization framework.

Additional outputs include:

```text
unnormalized_count.*
unnormalized_copy.*
normalized_16S.*
ppm.*
rpkm.*
tpm.*
```

## Database files

CARDOAP provides three main database files:

```text
database/
├── CARD_protein_homolog_full_6059.fasta
├── CARD_structure_full_6059.txt
└── CARD_OAP_full_6059_annotation.csv
```

`CARD_protein_homolog_full_6059.fasta` contains the 6,059 CARD protein homolog reference sequences.

`CARD_structure_full_6059.txt` provides the annotation structure used by ARGs-OAP.

`CARD_OAP_full_6059_annotation.csv` contains detailed CARD annotation metadata.

## MAG ARG module: CARD RGI + CARDOAP nCell

This optional module quantifies ARG ORFs predicted from MAGs. It uses the **same samples'** CARDOAP Stage 1 `nCell` estimates, but its MAG ORF counts are a separate result from the read-level `normalized_cell.*` tables above. Do not add the two totals: many reads can support both analyses.

1. Predict MAG ORFs and retain paired nucleotide (`representative.fna`) and protein (`representative.faa`) FASTA files with identical representative ORF IDs. If ORFs are clustered, retain the representative-to-original-ORF/MAG table. The RGI call below annotates representatives; cluster members are not independently validated by RGI.
2. Run CARD RGI main on the representative proteins using a local CARD database, and keep only `Perfect` and `Strict` matches. Preserve efflux pump hits unless a separate, documented sensitivity analysis removes them.
3. Map each sample's original, cleaned paired reads to **one combined nucleotide reference** of the selected MAG ORFs. Count CIGAR `M`, `=`, and `X` bases of primary mapped read mates per ORF. Divide by ORF length to obtain mean read-base depth, then by that sample's CARDOAP Stage 1 `nCell` to obtain copies/cell:

```text
MAG_ORF_copies_per_cell[sample, ORF]
  = total_aligned_bp[sample, ORF] / ORF_length_bp[ORF] / nCell[sample]
```

This preserves the original [HumandnohumARGs](https://github.com/uby76/HumandnohumARGs) calculation: Bowtie2 default end-to-end alignment with `--seed 0`, excluding unmapped, secondary and supplementary alignments (SAM flag mask `2308`), with no MAPQ threshold. R1 and R2 contribute separately, so overlapping mate bases count twice. PCR duplicates are not explicitly removed. These are read-base depths, not fragment-corrected depths.

### Run the MAG module

Install `rgi`, `bowtie2`, and `samtools`; use the same CARD/RGI release for all MAGs. First run RGI main on the representative proteins (repeat or parallelize for chunks if needed):

```bash
rgi main \
  --input_sequence representative.faa \
  --output_file mag_rgi/CARD_RGI \
  --input_type protein \
  --alignment_tool DIAMOND \
  --local --clean --data wgs
```

Make the nucleotide reference. `--rgi-dir` finds all `CARD_RGI.txt` files below the given directory. If representative ORFs were clustered, supply their original-ORF/MAG mapping with `--membership-table`. The HumandnohumARGs workflow produces `representative_to_original_ORF_MAG.tsv.gz` with the required `final_representative_ORF`, `original_ORF`, and `MAG_ID` columns. `--orf-metadata` is optional and may contain `ORF_ID`, `Source_category`, and `representative_MAG_ID` columns.

```bash
python3 scripts/mag_orf_copies_per_cell.py prepare \
  --representative-fna representative.fna \
  --rgi-dir mag_rgi \
  --membership-table representative_to_original_ORF_MAG.tsv.gz \
  --output-dir mag_arg_reference

bowtie2-build mag_arg_reference/CARD_MAG_ORFs.fna mag_arg_reference/index
```

Use the `sample` and `nCell` columns from `output/metadata.txt` produced by CARDOAP Stage 1. Sample IDs must match exactly. Run once per sample; use a job array for large cohorts. The same cleaned reads should be used for both CARDOAP Stage 1 and MAG ORF mapping.

```bash
python3 scripts/mag_orf_copies_per_cell.py quantify \
  --reference-dir mag_arg_reference \
  --cell-metadata output/metadata.txt \
  --sample sample1 \
  --read-one clean_reads/sample1_R1.fastq.gz \
  --read-two clean_reads/sample1_R2.fastq.gz \
  --threads 8 \
  --output-dir mag_arg_results/sample1
```

`mag_arg_results/sample1/sample1.CARD_MAG_ORF_copies_per_cell.tsv` contains sample × representative ARG ORF aligned bases, estimated gene copies, copies/cell, and the MAG ID when it can be assigned to one MAG. `sample1.CARD_MAG_by_MAG.tsv` sums only representative ORFs whose cluster maps to one MAG. `sample1.CARD_MAG_summary.tsv` separates the assigned and unassigned copies/cell. `mag_arg_reference/CARD_MAG_ORF_annotation.tsv` records the RGI evidence; `CARD_MAG_ORF_membership.tsv` records every original ORF ID and MAG ID associated with each selected representative. To reuse an existing BAM mapped to **this exact ORF reference**, pass `--bam sample.bam` instead of the read arguments.

For HumGut/NonHumGut analysis, the membership table supplies `Source_category` from its `category` column. Alternatively, provide `--orf-metadata representative_ORF_metadata.tsv` with one row per representative ORF and columns `ORF_ID`, `Source_category`, and, if available, `representative_MAG_ID`. Quantify against the same combined reference, then split the resulting ORF table by `Source_category`. Do not map the source groups separately.

An RGI result belongs to the **representative ORF**, whereas `CARD_MAG_ORF_membership.tsv` traces that representative to original MAG ORF IDs. If one representative contains multiple original ORFs, its read-based abundance cannot be copied to each original ORF without double counting. When all members belong to one MAG, the abundance can be included in that MAG's summary, but no separate value is assigned to each original ORF. When members span multiple MAGs, `MAG_ID` is blank and the abundance remains in the unassigned-cluster total. A per-original-gene estimate requires a separate mapping analysis with a reference and uniquely assignable reads at that resolution.

These results are named MAG ORF copies/cell. They are not equivalent to CARDOAP's read-level CARD reference copies/cell because the reference sequences and assignment procedure differ. `nCell` is a cell-equivalent denominator, not a 16S read count. Resfams-only candidates, if assessed later, should be reported separately from RGI Perfect/Strict ORFs.

## Workflow

```text
Short reads
    ↓
DIAMOND BLASTX
    ↓
Candidate ARG-like reads
    ↓
NCBI BLASTX
    ↓
CARD ARG identification
    ↓
Reference-length correction
    ↓
KO30 single-copy marker normalization
    ↓
ARG copies per cell
```

## Notes

The read-level CARDOAP workflow does not modify the core ARGs-OAP algorithm. It replaces the default SARG reference database with CARD protein homolog references while retaining the original ARGs-OAP search, filtering, and normalization framework. The optional MAG module above uses its `nCell` output as the denominator for an independent ORF read-mapping analysis.

CARD references can be associated with multiple drug classes, AMR gene families, or resistance mechanisms. Reference-, ARO-, and ARG-level results are therefore recommended for primary quantitative analyses.

## References

CARD:  
https://card.mcmaster.ca/

ARGs-OAP:  
https://github.com/xinehc/args_oap
