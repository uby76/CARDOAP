#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(dirname "$SCRIPT_DIR")"

cd "$ROOT_DIR"

DB="database/CARD_protein_homolog_full_6059.fasta"
STRUCTURE="database/CARD_structure_full_6059.txt"
THREADS=8

echo "Building CARDOAP database..."

args_oap make_db \
    -i "$DB"

echo "Running ARGs-OAP Stage 1..."

args_oap stage_one \
    -i test \
    -o test_output \
    -f fa \
    -t "$THREADS" \
    --database "$DB"

echo "Running ARGs-OAP Stage 2..."

args_oap stage_two \
    -i test_output \
    -t "$THREADS" \
    --database "$DB" \
    --structure1 "$STRUCTURE"

echo "CARDOAP test completed successfully."
