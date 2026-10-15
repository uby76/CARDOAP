#!/usr/bin/env python3
"""Predict ORFs in all MAGs and cluster them into one CARD-RGI input.

Input manifest: TSV with MAG_ID and fasta columns. Output retains every
representative-to-original-ORF-to-MAG link for downstream quantification.
"""

import argparse
import csv
import gzip
import json
import re
import sqlite3
import subprocess
from pathlib import Path


CD_HIT_PARAMETERS = ["-c", "0.95", "-aS", "0.90", "-G", "1", "-g", "1",
                     "-r", "1", "-n", "10", "-d", "0", "-l", "10"]
CLUSTER_RECORD = re.compile(r"\d+\s+\d+nt, >(.+)\.\.\. (.+)$")


def fasta_records(path):
    header, parts = None, []
    with Path(path).open() as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts).upper()
                header, parts = line[1:], []
            elif line:
                if header is None:
                    raise ValueError(f"Sequence before FASTA header: {path}")
                parts.append(line)
    if header is not None:
        yield header, "".join(parts).upper()


def clusters(path):
    members = []
    with Path(path).open() as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith(">Cluster"):
                if members:
                    yield members
                members = []
            elif line:
                match = CLUSTER_RECORD.fullmatch(line)
                if not match:
                    raise ValueError(f"Invalid CD-HIT cluster line: {line}")
                members.append((match.group(1), match.group(2) == "*"))
    if members:
        yield members


def read_manifest(path):
    path = Path(path)
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or not {"MAG_ID", "fasta"}.issubset(rows[0]):
        raise ValueError("Manifest must contain MAG_ID and fasta columns")
    ids = set()
    for row in rows:
        mid = row["MAG_ID"]
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", mid) or mid in (".", "..") or mid in ids:
            raise ValueError(f"Invalid or duplicate MAG_ID: {mid}")
        ids.add(mid)
        fasta = Path(row["fasta"])
        if not fasta.is_absolute():
            fasta = path.parent / fasta
        if not fasta.is_file():
            raise FileNotFoundError(fasta)
        row["fasta"] = str(fasta.resolve())
    return rows


def predict(manifest, out, database):
    all_dna = out / "all_ORFs.fna"
    all_protein = out / "all_ORFs.faa"
    fields = ["ORF_ID", "MAG_ID", "original_prodigal_ORF", "original_contig_ID",
              "start", "end", "strand", "ORF_length_bp"]
    count = 0
    with all_dna.open("w") as dna_out, all_protein.open("w") as protein_out, (out / "original_ORF_metadata.tsv").open("w", newline="") as meta_out:
        writer = csv.DictWriter(meta_out, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in manifest:
            mid = row["MAG_ID"]
            mag_dir = out / "prodigal" / mid
            mag_dir.mkdir(parents=True, exist_ok=True)
            fna, faa = mag_dir / "genes.ffn", mag_dir / "proteins.faa"
            cmd = ["prodigal", "-i", row["fasta"], "-d", str(fna), "-a", str(faa),
                   "-o", str(mag_dir / "genes.gff"), "-f", "gff", "-p", "meta", "-q"]
            with (mag_dir / "prodigal.log").open("w") as log:
                subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
            proteins = {}
            for header, sequence in fasta_records(faa):
                oid = header.split()[0]
                if oid in proteins:
                    raise ValueError(f"Duplicate Prodigal protein ID in {mid}: {oid}")
                proteins[oid] = sequence.rstrip("*")
            found = set()
            for header, sequence in fasta_records(fna):
                parts = header.split(" # ")
                oid = parts[0].split()[0]
                if len(parts) < 4 or oid in found or oid not in proteins:
                    raise ValueError(f"Invalid or unmatched Prodigal ORF in {mid}: {oid}")
                found.add(oid)
                unique = f"{mid}__{oid}"
                start, end, strand = int(parts[1]), int(parts[2]), int(parts[3])
                if (not sequence or not proteins[oid] or start < 1 or end < start or
                        strand not in (-1, 1) or len(sequence) != end - start + 1 or
                        re.search(r"[^ACGTRYSWKMBDHVN]", sequence)):
                    raise ValueError(f"Invalid ORF sequence or coordinates: {unique}")
                database.execute("INSERT INTO orfs VALUES (?, ?)", (unique, mid))
                dna_out.write(f">{unique}\n{sequence}\n")
                protein_out.write(f">{unique}\n{proteins[oid]}\n")
                writer.writerow({"ORF_ID": unique, "MAG_ID": mid, "original_prodigal_ORF": oid,
                                 "original_contig_ID": oid.rsplit("_", 1)[0],
                                 "start": start, "end": end, "strand": strand,
                                 "ORF_length_bp": len(sequence)})
                count += 1
            if not found or found != set(proteins):
                raise ValueError(f"Empty or mismatched Prodigal nucleotide/protein IDs: {mid}")
    database.commit()
    return count


def cluster(out, database, original_count, threads, memory_mb):
    cmd = ["cd-hit-est", "-i", str(out / "all_ORFs.fna"), "-o", str(out / "clustered.fna"),
           *CD_HIT_PARAMETERS, "-T", str(threads), "-M", str(memory_mb)]
    with (out / "cd_hit_est.log").open("w") as log:
        subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
    database.execute("CREATE TABLE members (id TEXT PRIMARY KEY)")
    representatives = set()
    links = 0
    mapping = out / "representative_to_original_ORF_MAG.tsv.gz"
    with gzip.open(mapping, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["final_representative_ORF", "original_ORF", "MAG_ID"])
        for group in clusters(out / "clustered.fna.clstr"):
            reps = [oid for oid, representative in group if representative]
            if len(reps) != 1:
                raise ValueError("CD-HIT cluster must contain exactly one representative")
            rep = reps[0]
            if rep in representatives:
                raise ValueError(f"Repeated representative ORF: {rep}")
            representatives.add(rep)
            for oid, _ in group:
                database.execute("INSERT INTO members VALUES (?)", (oid,))
                record = database.execute("SELECT mag FROM orfs WHERE id = ?", (oid,)).fetchone()
                if record is None:
                    raise ValueError(f"CD-HIT returned unknown ORF: {oid}")
                writer.writerow([rep, oid, record[0]])
                links += 1
    if links != original_count or database.execute("SELECT COUNT(*) FROM members").fetchone()[0] != original_count:
        raise ValueError("Not every original ORF occurs exactly once in the cluster membership")
    for source, target in (("all_ORFs.fna", "representative.fna"), ("all_ORFs.faa", "representative.faa")):
        seen = set()
        with (out / target).open("w") as handle:
            for header, sequence in fasta_records(out / source):
                oid = header.split()[0]
                if oid in representatives:
                    handle.write(f">{oid}\n{sequence}\n")
                    seen.add(oid)
        if seen != representatives:
            raise ValueError(f"Missing representative sequences in {source}")
    database.commit()
    return len(representatives), links


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, help="TSV columns MAG_ID and fasta")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--memory-mb", type=int, default=0, help="CD-HIT-EST memory limit; 0 means unlimited")
    args = parser.parse_args()
    if args.threads < 1 or args.memory_mb < 0:
        parser.error("--threads must be positive and --memory-mb nonnegative")
    manifest = read_manifest(args.manifest)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    database = sqlite3.connect(out / ".orf_membership_check.sqlite")
    try:
        database.execute("DROP TABLE IF EXISTS orfs")
        database.execute("DROP TABLE IF EXISTS members")
        database.execute("CREATE TABLE orfs (id TEXT PRIMARY KEY, mag TEXT NOT NULL)")
        originals = predict(manifest, out, database)
        reps, links = cluster(out, database, originals, args.threads, args.memory_mb)
    finally:
        database.close()
        (out / ".orf_membership_check.sqlite").unlink(missing_ok=True)
    (out / "prepare_summary.json").write_text(json.dumps({"MAGs": len(manifest), "original_ORFs": originals,
        "representative_ORFs": reps, "membership_links": links,
        "prodigal_mode": "meta", "cd_hit_est_parameters": CD_HIT_PARAMETERS}, indent=2) + "\n")
    print(f"Prepared {originals} ORFs from {len(manifest)} MAGs; {reps} representatives")


if __name__ == "__main__":
    main()
