```markdown
# CARDOAP

CARDOAP integrates the CARD protein homolog database with the ARGs-OAP framework for short-read ARG profiling and cell-normalized quantification.

The current database contains **6,059 CARD protein homolog reference sequences**. CARDOAP retains the original two-stage ARGs-OAP workflow and KO30-based cell normalization, while replacing the default SARG database with CARD.

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

CARDOAP does not modify the core ARGs-OAP algorithm. It replaces the default SARG reference database with CARD protein homolog references while retaining the original ARGs-OAP search, filtering, and normalization framework.

CARD references can be associated with multiple drug classes, AMR gene families, or resistance mechanisms. Reference-, ARO-, and ARG-level results are therefore recommended for primary quantitative analyses.

## References

CARD:  
https://card.mcmaster.ca/

ARGs-OAP:  
https://github.com/xinehc/args_oap
