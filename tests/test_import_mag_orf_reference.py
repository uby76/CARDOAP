import csv
import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/import_mag_orf_reference.py"


class ImportMagOrfReferenceTest(unittest.TestCase):
    def test_import_preserves_representatives_and_original_mag_ids(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            upstream = root / "reference"
            upstream.mkdir()
            for category, records in {
                "Human-specific": [("MAG_A__orf1", "ATGAAA", "MK")],
                "NonHuman-specific": [],
                "Shared": [("MAG_B__orf2", "ATGCCC", "MP")],
            }.items():
                (upstream / f"{category}.fna").write_text("".join(f">{oid}\n{dna}\n" for oid, dna, _ in records))
                (upstream / f"{category}.faa").write_text("".join(f">{oid}\n{aa}\n" for oid, _, aa in records))
            with gzip.open(upstream / "representative_to_original_ORF_MAG.tsv.gz", "wt") as handle:
                handle.write("final_representative_ORF\tcategory\toriginal_ORF\tMAG_ID\n"
                             "MAG_A__orf1\tHuman-specific\tMAG_A__orf1\tMAG_A\n"
                             "MAG_B__orf2\tShared\tMAG_B__orf2\tMAG_B\n"
                             "MAG_B__orf2\tShared\tMAG_C__orf7\tMAG_C\n")
            out = root / "imported"
            subprocess.run([sys.executable, str(SCRIPT), "--upstream-reference-dir", str(upstream),
                            "--output-dir", str(out)], check=True, capture_output=True)
            self.assertIn(">MAG_A__orf1", (out / "representative.fna").read_text())
            self.assertIn(">MAG_B__orf2", (out / "representative.faa").read_text())
            with gzip.open(out / "representative_to_original_ORF_MAG.tsv.gz", "rt") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([(r["original_ORF"], r["MAG_ID"]) for r in rows],
                             [("MAG_A__orf1", "MAG_A"), ("MAG_B__orf2", "MAG_B"), ("MAG_C__orf7", "MAG_C")])
            report = json.loads((out / "import_summary.json").read_text())
            self.assertEqual((report["representative_ORFs"], report["original_ORFs"]), (2, 3))


if __name__ == "__main__":
    unittest.main()
