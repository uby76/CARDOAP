#!/usr/bin/env python3
"""Quantify CARD RGI Perfect/Strict MAG representative ORFs per cell.

The aligned-base definition reproduces the HumandnohumARGs MAG workflow:
primary mapped reads, CIGAR M/= /X bases, and mates counted separately.
"""

import argparse
import csv
import gzip
import math
import re
import subprocess
from collections import Counter
from pathlib import Path


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
EXCLUDED_FLAGS = 4 | 256 | 2048


def table_rows(path):
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t" if ".tsv" in path.suffixes or ".txt" in path.suffixes else ",")
        if not reader.fieldnames:
            raise ValueError(f"Empty table: {path}")
        yield from reader


def read_table(path):
    return list(table_rows(path))


def write_table(path, fields, rows):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def fasta_records(path):
    name, parts = None, []
    with Path(path).open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts).upper()
                name, parts = line[1:].split()[0], []
            elif name is None:
                raise ValueError(f"Sequence before FASTA header: {path}")
            else:
                parts.append(line)
    if name is not None:
        yield name, "".join(parts).upper()


def prepare(args):
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    tables = [Path(p) for p in args.rgi_table]
    if args.rgi_dir:
        tables.extend(sorted(Path(args.rgi_dir).rglob("CARD_RGI.txt")))
    if not tables:
        raise ValueError("Provide --rgi-table or --rgi-dir")
    hits = {}
    for table in tables:
        for row in read_table(table):
            oid = row["ORF_ID"].split()[0]
            if row["Cut_Off"] not in ("Perfect", "Strict"):
                continue
            if oid in hits:
                previous = hits[oid]
                rank = lambda r: (r["Cut_Off"] == "Perfect", float(r.get("Best_Hit_Bitscore") or 0))
                if rank(row) < rank(previous):
                    continue
                if rank(row) == rank(previous) and row.get("ARO") != previous.get("ARO"):
                    raise ValueError(f"Ambiguous equal-rank RGI hits for representative ORF: {oid}")
            hits[oid] = row
    if not hits:
        raise ValueError("No CARD Perfect/Strict representative ORFs")
    metadata = {}
    if args.orf_metadata:
        for row in read_table(args.orf_metadata):
            oid = row["ORF_ID"]
            if oid in metadata:
                raise ValueError(f"Duplicate metadata ORF: {oid}")
            metadata[oid] = row
    membership = {}
    member_categories = {}
    member_rows = []
    if args.membership_table:
        original_seen = set()
        for row in table_rows(args.membership_table):
            rep = row["final_representative_ORF"]
            if rep not in hits:
                continue
            original, mag = row["original_ORF"], row["MAG_ID"]
            if not original or not mag or original in original_seen:
                raise ValueError(f"Invalid or duplicate original ORF in membership table: {original}")
            original_seen.add(original)
            membership.setdefault(rep, []).append((original, mag))
            member_categories.setdefault(rep, set()).add(row.get("category", ""))
            member_rows.append({"representative_ORF_ID": rep, "original_ORF_ID": original, "MAG_ID": mag,
                                "Source_category": row.get("category", ""),
                                "original_contig_ID": row.get("original_contig_ID", ""),
                                "start": row.get("start", ""), "end": row.get("end", ""),
                                "strand": row.get("strand", "")})
        if set(membership) != set(hits):
            raise ValueError(f"Membership missing for {len(set(hits) - set(membership))} selected representative ORFs")
        if any(len(categories) != 1 for categories in member_categories.values()):
            raise ValueError("Conflicting source categories for a representative ORF")
    seen, annotations = set(), []
    with (out / "CARD_MAG_ORFs.fna").open("w") as handle:
        for oid, seq in fasta_records(args.representative_fna):
            if oid in seen:
                raise ValueError(f"Duplicate FASTA ORF: {oid}")
            seen.add(oid)
            if oid not in hits:
                continue
            if not seq or re.search(r"[^ACGTRYSWKMBDHVN]", seq):
                raise ValueError(f"Invalid nucleotide sequence for {oid}")
            hit = hits[oid]
            source = metadata.get(oid, {})
            members = membership.get(oid, [])
            singleton = len(members) == 1
            mags = {mag for _, mag in members}
            source_category = source.get("Source_category", source.get("category", ""))
            membership_category = next(iter(member_categories[oid])) if members else ""
            if source_category and membership_category and source_category != membership_category:
                raise ValueError(f"Conflicting source category for {oid}")
            representative_mag = source.get("representative_MAG_ID", "")
            if not representative_mag:
                representative_mag = next((mag for original, mag in members if original == oid), "")
            annotations.append({
                "ORF_ID": oid,
                "ORF_length_bp": len(seq),
                "Cut_Off": hit["Cut_Off"],
                "Best_Hit_ARO": hit.get("Best_Hit_ARO", ""),
                "ARO": hit.get("ARO", ""),
                "Source_category": source_category or membership_category,
                "representative_MAG_ID": representative_mag,
                "MAG_ID": next(iter(mags)) if len(mags) == 1 else "",
                "original_ORF_ID": members[0][0] if singleton else "",
                "original_ORF_count": len(members),
                "member_MAG_count": len(mags),
                "MAG_assignment": ("unique_original_ORF" if singleton else
                                   "unique_MAG_cluster" if len(mags) == 1 else
                                   "multiple_MAG_cluster" if members else "membership_not_supplied"),
            })
            handle.write(f">{oid}\n")
            for pos in range(0, len(seq), 80):
                handle.write(seq[pos:pos + 80] + "\n")
    missing = set(hits) - seen
    if missing:
        raise ValueError(f"RGI hit ORFs absent from nucleotide FASTA: {len(missing)} (e.g. {min(missing)})")
    write_table(out / "CARD_MAG_ORF_annotation.tsv", list(annotations[0]), annotations)
    if args.membership_table:
        write_table(out / "CARD_MAG_ORF_membership.tsv", list(member_rows[0]), member_rows)
    print(f"Prepared {len(annotations)} CARD Perfect/Strict representative ORFs in {out}")


def cigar_aligned_bp(cigar):
    terms = CIGAR_RE.findall(cigar)
    if not terms or "".join(n + op for n, op in terms) != cigar:
        raise ValueError(f"Invalid CIGAR: {cigar}")
    return sum(int(n) for n, op in terms if op in "M=X")


def ncell_for_sample(path, sample):
    rows = read_table(path)
    matches = [r for r in rows if r.get("sample") == sample]
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one nCell row for {sample}; found {len(matches)}")
    value = float(matches[0]["nCell"])
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"Invalid nCell for {sample}: {value}")
    return value


def bam_aligned_bp(bam, lengths):
    subprocess.run(["samtools", "quickcheck", str(bam)], check=True)
    header = subprocess.run(["samtools", "view", "-H", str(bam)], check=True, capture_output=True, text=True).stdout
    observed = {}
    for line in header.splitlines():
        if line.startswith("@SQ\t"):
            tags = dict(item.split(":", 1) for item in line.split("\t")[1:])
            observed[tags["SN"]] = int(tags["LN"])
    if observed != lengths:
        raise ValueError("BAM reference names/lengths differ from prepared MAG ORFs")
    bases, reads = Counter(), Counter()
    process = subprocess.Popen(["samtools", "view", "-F", str(EXCLUDED_FLAGS), str(bam)], stdout=subprocess.PIPE, text=True)
    try:
        for line in process.stdout:
            fields = line.rstrip("\n").split("\t")
            oid = fields[2]
            bases[oid] += cigar_aligned_bp(fields[5])
            reads[oid] += 1
    finally:
        process.stdout.close()
    if process.wait():
        raise RuntimeError("samtools view failed")
    return bases, reads


def quantify(args):
    ref = Path(args.reference_dir)
    annotations = read_table(ref / "CARD_MAG_ORF_annotation.tsv")
    lengths = {r["ORF_ID"]: int(r["ORF_length_bp"]) for r in annotations}
    if not lengths or len(lengths) != len(annotations) or any(x <= 0 for x in lengths.values()):
        raise ValueError("Invalid MAG ORF annotation table")
    cell = ncell_for_sample(args.cell_metadata, args.sample)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    bam = out / f"{args.sample}.primary_mapped.bam"
    if not args.bam:
        if not all((ref / f"index.{suffix}").exists() for suffix in ("1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2")) and not all((ref / f"index.{suffix}").exists() for suffix in ("1.bt2l", "2.bt2l", "3.bt2l", "4.bt2l", "rev.1.bt2l", "rev.2.bt2l")):
            raise FileNotFoundError(f"Bowtie2 index absent: {ref}/index.*; run bowtie2-build first")
        command = ["bowtie2", "--seed", "0", "-x", str(ref / "index"), "-p", str(args.threads), "-1", args.read_one, "-2", args.read_two]
        with (out / f"{args.sample}.bowtie2.log").open("w") as log, (out / f"{args.sample}.samtools.log").open("w") as err:
            mapper = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=log)
            converter = subprocess.Popen(["samtools", "view", "-b", "-F", str(EXCLUDED_FLAGS), "-o", str(bam), "-"], stdin=mapper.stdout, stderr=err)
            mapper.stdout.close()
            convert_code, map_code = converter.wait(), mapper.wait()
        if map_code or convert_code:
            raise RuntimeError(f"Bowtie2/samtools failed: {map_code}/{convert_code}")
    else:
        bam = Path(args.bam)
    bases, reads = bam_aligned_bp(bam, lengths)
    rows = []
    for ann in annotations:
        oid, length = ann["ORF_ID"], lengths[ann["ORF_ID"]]
        rows.append({**ann, "sample": args.sample, "nCell": format(cell, ".17g"),
                     "mapped_read_mates": reads[oid], "total_aligned_bp": bases[oid],
                     "estimated_gene_copies": format(bases[oid] / length, ".17g"),
                     "copies_per_cell": format(bases[oid] / length / cell, ".17g")})
    fields = ["sample", "nCell"] + list(annotations[0]) + ["mapped_read_mates", "total_aligned_bp", "estimated_gene_copies", "copies_per_cell"]
    write_table(out / f"{args.sample}.CARD_MAG_ORF_copies_per_cell.tsv", fields, rows)
    total = sum(float(row["copies_per_cell"]) for row in rows)
    by_mag = {}
    for row in rows:
        mag = row["MAG_ID"]
        if mag:
            count, amount = by_mag.get(mag, (0, 0.0))
            by_mag[mag] = (count + 1, amount + float(row["copies_per_cell"]))
    write_table(out / f"{args.sample}.CARD_MAG_by_MAG.tsv",
                ["sample", "MAG_ID", "representative_ARG_ORFs", "ARG_copies_per_cell"],
                ({"sample": args.sample, "MAG_ID": mag, "representative_ARG_ORFs": count,
                  "ARG_copies_per_cell": format(amount, ".17g")}
                 for mag, (count, amount) in sorted(by_mag.items())))
    assigned = sum(amount for _, amount in by_mag.values())
    unassigned = sum(float(row["copies_per_cell"]) for row in rows if not row["MAG_ID"])
    write_table(out / f"{args.sample}.CARD_MAG_summary.tsv",
                ["sample", "nCell", "CARD_MAG_ORFs", "CARD_MAG_total_copies_per_cell",
                 "assigned_to_one_MAG_copies_per_cell", "unassigned_cluster_copies_per_cell"],
                [{"sample": args.sample, "nCell": format(cell, ".17g"), "CARD_MAG_ORFs": len(rows),
                  "CARD_MAG_total_copies_per_cell": format(total, ".17g"),
                  "assigned_to_one_MAG_copies_per_cell": format(assigned, ".17g"),
                  "unassigned_cluster_copies_per_cell": format(unassigned, ".17g")}])
    print(f"{args.sample}: {len(rows)} ORFs; total MAG ARG copies/cell = {total:.9g}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prep = commands.add_parser("prepare", help="Select RGI Perfect/Strict ORFs and make nucleotide reference")
    prep.add_argument("--representative-fna", required=True, help="Nucleotide FASTA matching the RGI representative proteins")
    prep.add_argument("--rgi-table", action="append", default=[], help="RGI CARD_RGI.txt; may be repeated")
    prep.add_argument("--rgi-dir", help="Directory searched recursively for CARD_RGI.txt")
    prep.add_argument("--orf-metadata", help="Optional representative ORF annotation CSV/TSV")
    prep.add_argument("--membership-table", help="Representative-to-original-ORF/MAG TSV or TSV.GZ from MAG workflow")
    prep.add_argument("--output-dir", required=True)
    prep.set_defaults(func=prepare)
    quant = commands.add_parser("quantify", help="Map one sample and divide ORF mean depth by its CARDOAP nCell")
    quant.add_argument("--reference-dir", required=True)
    quant.add_argument("--cell-metadata", required=True, help="CARDOAP stage_one/metadata.txt with sample and nCell")
    quant.add_argument("--sample", required=True)
    source = quant.add_mutually_exclusive_group(required=True)
    source.add_argument("--bam", help="Existing Bowtie2 BAM against the prepared reference")
    source.add_argument("--read-one", help="Read 1 FASTQ/FASTA, optionally gzip compressed")
    quant.add_argument("--read-two", help="Read 2 FASTQ/FASTA")
    quant.add_argument("--threads", type=int, default=8)
    quant.add_argument("--output-dir", required=True)
    quant.set_defaults(func=quantify)
    args = parser.parse_args()
    if args.command == "quantify" and bool(args.read_one) != bool(args.read_two):
        parser.error("--read-one and --read-two must be supplied together")
    args.func(args)


if __name__ == "__main__":
    main()
