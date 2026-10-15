#!/usr/bin/env python3
"""Import the exact clustered MAG ORF reference from HumandnohumARGs.

Upstream: Prodigal -p meta, globally unique MAG__ORF IDs, within-group
CD-HIT-EST 95/90 clustering, bidirectional Human/NonHuman comparison,
and shared-ORF reclustering. This importer does not recluster sequences.
"""

import argparse
import csv
import gzip
import json
import sqlite3
from pathlib import Path

from mag_orf_copies_per_cell import fasta_records


CATEGORIES = ("Human-specific", "NonHuman-specific", "Shared")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upstream-reference-dir", required=True,
                        help="HumandnohumARGs 07_HumGut_ORF_CARD/reference directory")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    source, out = Path(args.upstream_reference_dir), Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    members = source / "representative_to_original_ORF_MAG.tsv.gz"
    if not members.is_file():
        raise FileNotFoundError(members)
    selected = {}
    with (out / "representative.fna").open("w") as nucleotide, (out / "representative.faa").open("w") as protein:
        for category in CATEGORIES:
            fna, faa = source / f"{category}.fna", source / f"{category}.faa"
            dna = list(fasta_records(fna))
            amino = list(fasta_records(faa))
            if [oid for oid, _ in dna] != [oid for oid, _ in amino]:
                raise ValueError(f"Nucleotide and protein representative IDs differ: {category}")
            for (oid, seq), (_, aa) in zip(dna, amino):
                if oid in selected or not seq or not aa:
                    raise ValueError(f"Duplicate or empty representative ORF: {oid}")
                selected[oid] = category
                nucleotide.write(f">{oid}\n{seq}\n")
                protein.write(f">{oid}\n{aa}\n")
    if not selected:
        raise ValueError("No representative MAG ORFs")

    database = sqlite3.connect(out / ".membership_check.sqlite")
    database.execute("CREATE TABLE originals (original TEXT PRIMARY KEY, representative TEXT, mag TEXT)")
    counts = {oid: 0 for oid in selected}
    copied = out / "representative_to_original_ORF_MAG.tsv.gz"
    try:
        with gzip.open(members, "rt", newline="") as source_handle, gzip.open(copied, "wt", newline="") as dest_handle:
            reader = csv.DictReader(source_handle, delimiter="\t")
            required = {"final_representative_ORF", "original_ORF", "MAG_ID", "category"}
            if not reader.fieldnames or not required.issubset(reader.fieldnames):
                raise ValueError(f"Membership table lacks {required}")
            writer = csv.DictWriter(dest_handle, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in reader:
                rep, original, mag = row["final_representative_ORF"], row["original_ORF"], row["MAG_ID"]
                if rep not in selected or row["category"] != selected[rep] or not original or not mag:
                    raise ValueError(f"Invalid representative/member/category mapping: {rep} / {original}")
                database.execute("INSERT INTO originals VALUES (?, ?, ?)", (original, rep, mag))
                counts[rep] += 1
                writer.writerow(row)
        if any(value == 0 for value in counts.values()):
            raise ValueError("Some representatives have no original ORF/MAG members")
        database.commit()
        originals = database.execute("SELECT COUNT(*) FROM originals").fetchone()[0]
    finally:
        database.close()
        (out / ".membership_check.sqlite").unlink(missing_ok=True)
    (out / "import_summary.json").write_text(json.dumps({"representative_ORFs": len(selected),
        "original_ORFs": originals, "by_category": {c: sum(v == c for v in selected.values()) for c in CATEGORIES},
        "source_reference_dir": str(source.resolve())}, indent=2) + "\n")
    print(f"Imported {len(selected)} representative ORFs and {originals} original ORF/MAG links")


if __name__ == "__main__":
    main()
