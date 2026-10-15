import csv
import gzip
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/mag_orf_copies_per_cell.py"


class MagOrfWorkflowTest(unittest.TestCase):
    def test_prepare_and_existing_bam_quantification(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "representatives.fna").write_text(">orfA\n" + "A" * 100 + "\n>orfB\n" + "C" * 80 + "\n>orfC\n" + "G" * 100 + "\n")
            (root / "CARD_RGI.txt").write_text(
                "ORF_ID\tCut_Off\tBest_Hit_ARO\tARO\tBest_Hit_Bitscore\n"
                "orfA\tStrict\tARG_A\t3000001\t90\n"
                "orfA\tPerfect\tARG_A\t3000001\t100\n"
                "orfB\tLoose\tARG_B\t3000002\t80\n"
                "orfC\tStrict\tARG_C\t3000003\t95\n"
            )
            with gzip.open(root / "members.tsv.gz", "wt") as handle:
                handle.write("final_representative_ORF\toriginal_ORF\tMAG_ID\n"
                             "orfA\torfA\tMAG_A\n"
                             "orfC\torfC\tMAG_C\n"
                             "orfC\torfD\tMAG_D\n")
            (root / "metadata.txt").write_text("sample\tnRead\tn16S\tnCell\nS1\t100\t5\t2\nS2\t100\t5\t4\n")
            ref = root / "ref"
            subprocess.run([sys.executable, str(SCRIPT), "prepare", "--representative-fna", str(root / "representatives.fna"),
                            "--rgi-table", str(root / "CARD_RGI.txt"), "--membership-table", str(root / "members.tsv.gz"),
                            "--output-dir", str(ref)], check=True, capture_output=True)
            with (ref / "CARD_MAG_ORF_annotation.tsv").open() as handle:
                annotation = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([(r["ORF_ID"], r["Cut_Off"]) for r in annotation], [("orfA", "Perfect"), ("orfC", "Strict")])
            self.assertEqual((annotation[0]["MAG_ID"], annotation[0]["original_ORF_ID"]), ("MAG_A", "orfA"))
            self.assertEqual(annotation[1]["MAG_assignment"], "multiple_MAG_cluster")
            self.assertEqual(annotation[1]["MAG_ID"], "")
            with (ref / "CARD_MAG_ORF_membership.tsv").open() as handle:
                self.assertEqual(len(list(csv.DictReader(handle, delimiter="\t"))), 3)

            # Fake samtools keeps this test independent of installed bioinformatics tools.
            # Two 20M mates overlap by 10 bp. The historical read-base rule counts 40 bp.
            bam = root / "mapped.bam"
            bam.write_text("@SQ\tSN:orfA\tLN:100\n@SQ\tSN:orfC\tLN:100\n"
                           "pair\t65\torfA\t1\t0\t20M\t*\t0\t0\t*\t*\n"
                           "pair\t129\torfA\t11\t0\t20M\t*\t0\t0\t*\t*\n"
                           "shared\t0\torfC\t1\t0\t10M\t*\t0\t0\t*\t*\n"
                           "secondary\t256\torfA\t1\t0\t30M\t*\t0\t0\t*\t*\n")
            bindir = root / "bin"
            bindir.mkdir()
            fake = bindir / "samtools"
            fake.write_text("#!/usr/bin/env python3\nimport sys\nfrom pathlib import Path\n"
                            "a=sys.argv[1:]\n"
                            "if a[0]=='quickcheck': sys.exit(0)\n"
                            "if a[0]=='view':\n"
                            "  lines=Path(a[-1]).read_text().splitlines()\n"
                            "  if '-H' in a: print('\\n'.join(x for x in lines if x.startswith('@'))); sys.exit(0)\n"
                            "  mask=int(a[a.index('-F')+1])\n"
                            "  print('\\n'.join(x for x in lines if not x.startswith('@') and not int(x.split('\\t')[1])&mask))\n"
                            "  sys.exit(0)\n"
                            "sys.exit(1)\n")
            fake.chmod(0o755)
            env = dict(os.environ, PATH=str(bindir) + os.pathsep + os.environ.get("PATH", ""))
            for sample, expected in (("S1", 0.2), ("S2", 0.1)):
                out = root / sample
                subprocess.run([sys.executable, str(SCRIPT), "quantify", "--reference-dir", str(ref),
                                "--cell-metadata", str(root / "metadata.txt"), "--sample", sample,
                                "--bam", str(bam), "--output-dir", str(out)], check=True, capture_output=True, env=env)
                with (out / f"{sample}.CARD_MAG_ORF_copies_per_cell.tsv").open() as handle:
                    rows = list(csv.DictReader(handle, delimiter="\t"))
                row = rows[0]
                self.assertEqual(int(row["total_aligned_bp"]), 40)
                self.assertEqual(int(row["mapped_read_mates"]), 2)
                self.assertAlmostEqual(float(row["estimated_gene_copies"]), 0.4)
                self.assertAlmostEqual(float(row["copies_per_cell"]), expected)
                self.assertEqual(rows[1]["MAG_ID"], "")
                with (out / f"{sample}.CARD_MAG_by_MAG.tsv").open() as handle:
                    mag_rows = list(csv.DictReader(handle, delimiter="\t"))
                self.assertEqual(len(mag_rows), 1)
                self.assertEqual(mag_rows[0]["MAG_ID"], "MAG_A")
                self.assertAlmostEqual(float(mag_rows[0]["ARG_copies_per_cell"]), expected)
                with (out / f"{sample}.CARD_MAG_summary.tsv").open() as handle:
                    summary = next(csv.DictReader(handle, delimiter="\t"))
                self.assertAlmostEqual(float(summary["unassigned_cluster_copies_per_cell"]), 0.1 / float(summary["nCell"]))


if __name__ == "__main__":
    unittest.main()
