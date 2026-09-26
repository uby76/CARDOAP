# CARDOAP

CARDOAP is a CARD-based antimicrobial resistance gene (ARG) database adapted for short-read ARG profiling and cell-normalized quantification using the ARGs-OAP framework.

CARDOAP replaces the default SARG database in ARGs-OAP with protein homolog reference sequences from the Comprehensive Antibiotic Resistance Database (CARD), while retaining the original ARGs-OAP two-stage search, filtering, and normalization framework.

## Database

The current CARDOAP database contains **6,059 CARD protein homolog model reference sequences**.

The database retains the complete CARD protein homolog reference set used for this release.

Each reference is linked to CARD annotation information including:

- ARO accession
- ARG name
- AMR gene family
- drug class
- resistance mechanism

The database files are:

```text
database/CARD_protein_homolog_full_6059.fasta
database/CARD_structure_full_6059.txt
database/CARD_OAP_full_6059_annotation.csv
