import csv
import gzip
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/prepare_mag_orfs.py"
QUANTIFY = Path(__file__).resolve().parents[1] / "scripts/mag_orf_copies_per_cell.py"


class PrepareMagOrfsTest(unittest.TestCase):
    def test_all_mags_are_clustered_with_original_ids_retained(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "MAG_A.fna").write_text(">contig\nATGAAA\n")
            (root / "MAG_B.fna").write_text(">contig\nATGAAA\n")
            (root / "mags.tsv").write_text("MAG_ID\tfasta\nMAG_A\tMAG_A.fna\nMAG_B\tMAG_B.fna\n")
            bin_dir = root / "bin"
            bin_dir.mkdir()
            prodigal = bin_dir / "prodigal"
            prodigal.write_text(
                "#!/usr/bin/env python3\n"
                "import sys\nfrom pathlib import Path\n"
                "a=sys.argv[1:]\n"
                "def value(x): return a[a.index(x)+1]\n"
                "Path(value('-d')).write_text('>contig_1 # 1 # 6 # 1 # ID=1\\nATGAAA\\n')\n"
                "Path(value('-a')).write_text('>contig_1 # 1 # 6 # 1 # ID=1\\nMK\\n')\n"
                "Path(value('-o')).write_text('##gff-version 3\\n')\n")
            prodigal.chmod(0o755)
            cdhit = bin_dir / "cd-hit-est"
            cdhit.write_text(
                "#!/usr/bin/env python3\n"
                "import sys\nfrom pathlib import Path\n"
                "a=sys.argv[1:]\n"
                "def value(x): return a[a.index(x)+1]\n"
                "out=Path(value('-o'))\n"
                "out.write_text('>MAG_A__contig_1\\nATGAAA\\n')\n"
                "Path(str(out)+'.clstr').write_text('>Cluster 0\\n0 6nt, >MAG_A__contig_1... *\\n1 6nt, >MAG_B__contig_1... at 100.00%\\n')\n")
            cdhit.chmod(0o755)
            env = dict(os.environ, PATH=str(bin_dir) + os.pathsep + os.environ.get("PATH", ""))
            out = root / "mag_orfs"
            subprocess.run([sys.executable, str(SCRIPT), "--manifest", str(root / "mags.tsv"),
                            "--output-dir", str(out), "--threads", "2", "--memory-mb", "128"],
                           check=True, capture_output=True, env=env)
            with gzip.open(out / "representative_to_original_ORF_MAG.tsv.gz", "rt") as handle:
                links = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([(r["original_ORF"], r["MAG_ID"]) for r in links],
                             [("MAG_A__contig_1", "MAG_A"), ("MAG_B__contig_1", "MAG_B")])
            self.assertEqual((out / "representative.fna").read_text(), ">MAG_A__contig_1\nATGAAA\n")
            self.assertEqual((out / "representative.faa").read_text(), ">MAG_A__contig_1\nMK\n")
            report = json.loads((out / "prepare_summary.json").read_text())
            self.assertEqual((report["MAGs"], report["original_ORFs"], report["representative_ORFs"]), (2, 2, 1))
            (root / "CARD_RGI.txt").write_text("ORF_ID\tCut_Off\tBest_Hit_ARO\tARO\n"
                                                  "MAG_A__contig_1\tStrict\tARG_A\t3000001\n")
            selected = root / "selected"
            subprocess.run([sys.executable, str(QUANTIFY), "prepare", "--representative-fna",
                            str(out / "representative.fna"), "--rgi-table", str(root / "CARD_RGI.txt"),
                            "--membership-table", str(out / "representative_to_original_ORF_MAG.tsv.gz"),
                            "--output-dir", str(selected)], check=True, capture_output=True)
            with (selected / "CARD_MAG_ORF_annotation.tsv").open() as handle:
                annotation = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(annotation["MAG_assignment"], "multiple_MAG_cluster")
            self.assertEqual(annotation["MAG_ID"], "")


if __name__ == "__main__":
    unittest.main()
