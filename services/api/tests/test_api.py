import json
import io
import sys
import zipfile
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from app import Handler, demo_data, derive
from analysis import compare_genomes, jaccard_similarity, pangenome_statistics, similarity_matrix, upgma
from ingest import parse_fasta, parse_gff3, parse_panaroo_csv
from io import StringIO
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "scripts"))
from ingest_real import parse_json_lines, safe_extract, select_records


class DemoDataTests(unittest.TestCase):
    def test_http_responses_include_security_headers(self):
        handler = object.__new__(Handler)
        headers = {}
        handler.send_response = lambda status: None
        handler.send_header = lambda name, value: headers.__setitem__(name, value)
        handler.end_headers = lambda: None
        handler.wfile = io.BytesIO()
        handler.respond(200, b"{}", "application/json")
        self.assertEqual(headers["X-Content-Type-Options"], "nosniff")
        self.assertEqual(headers["X-Frame-Options"], "DENY")
        self.assertEqual(headers["Referrer-Policy"], "no-referrer")
        self.assertIn("frame-ancestors 'none'", headers["Content-Security-Policy"])

    def test_fixture_is_deterministic(self):
        self.assertEqual(json.dumps(demo_data(), sort_keys=True), json.dumps(demo_data(), sort_keys=True))

    def test_derived_classes_and_copy_counts_are_consistent(self):
        result = derive(demo_data())
        self.assertEqual(sum(result["counts"].values()), len(result["families"]))
        self.assertEqual(result["counts"]["core"], 38)
        for strain in result["strains"]:
            genes = sum(f["copies"].get(strain["id"], 0) for f in result["families"])
            self.assertEqual(strain["genes"], genes)
            self.assertEqual(strain["families"], sum(strain["id"] in f["members"] for f in result["families"]))

    def test_prevalence_is_fraction_of_distinct_genomes(self):
        result = derive(demo_data())
        for family in result["families"]:
            self.assertAlmostEqual(family["prevalence"], len(family["members"]) / len(result["strains"]))

    def test_jaccard_similarity_uses_presence_not_copy_count(self):
        self.assertEqual(jaccard_similarity({"a", "b"}, {"b", "c"}), 1 / 3)
        self.assertEqual(jaccard_similarity(set(), set()), 1.0)

    def test_upgma_contains_each_genome_once_and_exports_newick(self):
        result = upgma(derive(demo_data()))
        self.assertIn("not a sequence-based", result["interpretation"])
        self.assertTrue(result["newick"].endswith(";"))
        leaf_names = []
        def visit(node):
            if node.get("name") is not None:
                leaf_names.append(node["name"])
            else:
                for child in node["children"]:
                    visit(child["node"])
        visit(result["root"])
        self.assertEqual(sorted(leaf_names), sorted(s["id"] for s in derive(demo_data())["strains"]))

    def test_statistics_are_deterministic_and_frequency_counts_balance(self):
        dataset = derive(demo_data())
        a = pangenome_statistics(dataset, permutations=12, seed=42)
        b = pangenome_statistics(dataset, permutations=12, seed=42)
        self.assertEqual(a, b)
        self.assertEqual(sum(row["families"] for row in a["family_frequency"]), len(dataset["families"]))
        self.assertEqual(a["accumulation"][-1]["pangenome_mean"], len(dataset["families"]))
        selected = pangenome_statistics(dataset, genome_ids=["LP-01", "LP-02"], permutations=8, seed=7)
        self.assertEqual(selected["selection_mode"], "selected genomes")
        self.assertEqual(selected["genomes"], ["LP-01", "LP-02"])
        self.assertEqual(selected["accumulation"][-1]["pangenome_mean"], len(set(f["id"] for f in dataset["families"] if set(f["members"]) & {"LP-01", "LP-02"})))

    def test_comparison_validates_selection_and_reports_shared_families(self):
        result = compare_genomes(derive(demo_data()), ["LP-01", "LP-02"])
        self.assertLessEqual(result["shared_family_count"], result["union_family_count"])
        self.assertIn("GF0001", result["family_ids"])
        with self.assertRaises(ValueError):
            compare_genomes(derive(demo_data()), ["LP-01"])

    def test_fasta_and_gff3_coordinate_normalization_and_bounds(self):
        sequences = parse_fasta(StringIO(">ctg one\nACGT\nACGN\n>ctg2\nTTTT\n"))
        genes = parse_gff3(StringIO("##gff-version 3\nctg\tNCBI\tCDS\t2\t5\t.\t-\t0\tID=cds-1;product=hypothetical%20protein\n"), {key: len(value) for key, value in sequences.items()})
        self.assertEqual(genes[0]["start"], 1)
        self.assertEqual(genes[0]["end"], 5)
        self.assertEqual(genes[0]["strand"], "-")
        self.assertEqual(genes[0]["attributes"]["product"], "hypothetical protein")
        with self.assertRaises(ValueError):
            parse_gff3(StringIO("ctg\tNCBI\tCDS\t1\t99\t.\t+\t0\tID=x\n"), {"ctg": 4})
        multipart = parse_gff3(StringIO("ctg\tNCBI\tCDS\t1\t2\t.\t+\t0\tID=cds-x;Parent=gene-x\nctg\tNCBI\tCDS\t3\t4\t.\t+\t1\tID=cds-x;Parent=gene-x\n"), {"ctg": 4})
        self.assertEqual(len(multipart), 2)
        self.assertEqual({gene["id"] for gene in multipart}, {"cds-x"})

    def test_panaroo_roary_csv_import_preserves_paralog_copy_counts(self):
        csv_text = "Gene,Non-unique Gene name,Annotation,No. isolates,No. sequences,Avg sequences per isolate,Genome Fragment,Order within Fragment,Accessory Fragment,Accessory Order within Fragment,QC,sample_A,sample_B\nGF1,,enzyme,2,3,1.5,,,,,,a1;a2,b1\nGF2,,unknown,1,1,1.0,,,,,,c1,\n"
        result = parse_panaroo_csv(StringIO(csv_text), {"sample_A": "GCF_A", "sample_B": "GCF_B"})
        self.assertEqual(result["families"][0]["copies"], {"GCF_A": 2, "GCF_B": 1})
        self.assertEqual(result["families"][1]["members"], ["GCF_A"])
        self.assertEqual(result["warnings"], [])

    def test_panaroo_import_rejects_missing_genome_columns(self):
        with self.assertRaises(ValueError):
            parse_panaroo_csv(StringIO("Gene,Annotation\nGF1,unknown\n"))

    def test_panaroo_18_matrix_without_roary_metadata_columns(self):
        result = parse_panaroo_csv(StringIO("Gene,Non-unique Gene name,Annotation,sample_A,sample_B\nGF1,,enzyme,a1;a2,\n"))
        self.assertEqual(result["families"][0]["members"], ["sample_A"])
        self.assertEqual(result["families"][0]["copies"], {"sample_A": 2})

    def test_ncbi_selection_ranks_quality_deduplicates_and_excludes_clinical_sources(self):
        def record(accession, level, strain, contigs, n50, source="fermented plant"):
            return {"accession": accession, "organism": {"organism_name": "Lactiplantibacillus plantarum", "infraspecific_names": {"strain": strain}}, "assembly_info": {"assembly_level": level, "total_sequence_length": "4000000", "number_of_contigs": str(contigs), "contig_n50": str(n50)}, "biosample": {"isolation_source": source}}
        reports = [record("GCF_000000001.1", "Scaffold", "strain A", 10, 100), record("GCF_000000002.1", "Complete Genome", "strain A", 1, 4000000), record("GCF_000000003.1", "Complete Genome", "strain B", 1, 4000000, "hospital patient"), record("GCF_000000004.1", "Chromosome", "strain C", 2, 2000000)]
        selected, excluded = select_records(reports, 3)
        self.assertEqual([row["accession"] for row in selected], ["GCF_000000002.1", "GCF_000000004.1"])
        self.assertEqual(len(excluded), 2)
        self.assertEqual(parse_json_lines('{"reports":[{"accession":"GCF_000000001.1"}]}')[0]["accession"], "GCF_000000001.1")

    def test_ncbi_package_extraction_rejects_path_traversal(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            package = base / "unsafe.zip"
            with zipfile.ZipFile(package, "w") as archive:
                archive.writestr("../escape.txt", "bad")
            with self.assertRaises(RuntimeError):
                safe_extract(package, base / "out")

    def test_api_get_routes_and_exports(self):
        from app import load_data
        routes = [
            ("/api/health", 200, "application/json"),
            ("/api/datasets", 200, "application/json"),
            ("/api/overview?dataset=demo-lp-v1", 200, "application/json"),
            ("/api/matrix?dataset=demo-lp-v1", 200, "application/json"),
            ("/api/strains?dataset=demo-lp-v1", 200, "application/json"),
            ("/api/families?dataset=demo-lp-v1", 200, "application/json"),
            ("/api/similarity?dataset=demo-lp-v1", 200, "application/json"),
            ("/api/tree?dataset=demo-lp-v1", 200, "application/json"),
            ("/api/statistics?dataset=demo-lp-v1&genomes=", 200, "application/json"),
            ("/api/statistics?dataset=demo-lp-v1&genomes=LP-01,LP-02", 200, "application/json"),
            ("/api/export/families.tsv?dataset=demo-lp-v1", 200, "text/tab-separated-values"),
            ("/api/export/families.csv?dataset=demo-lp-v1", 200, "text/csv"),
            ("/api/export/tree.newick?dataset=demo-lp-v1", 200, "text/plain"),
            ("/api/overview?dataset=nonexistent-dataset", 404, "application/json"),
            ("/", 200, "text/html"),
            ("/styles.css", 200, "text/css"),
            ("/app.js", 200, "application/javascript"),
        ]
        for path, expected_status, expected_mime in routes:
            handler = object.__new__(Handler)
            handler.path = path
            headers = {}
            status_holder = []
            handler.send_response = lambda s: status_holder.append(s)
            handler.send_header = lambda name, value: headers.__setitem__(name, value)
            handler.end_headers = lambda: None
            handler.wfile = io.BytesIO()
            handler.do_GET()
            self.assertEqual(status_holder[0], expected_status, f"Route {path} returned status {status_holder[0]}, expected {expected_status}")
            self.assertIn(expected_mime, headers["Content-Type"])

    def test_api_post_compare(self):
        handler = object.__new__(Handler)
        handler.path = "/api/compare"
        payload = json.dumps({"genome_ids": ["LP-01", "LP-02"], "dataset_id": "demo-lp-v1"}).encode()
        headers = {}
        status_holder = []
        handler.headers = {"Content-Length": str(len(payload))}
        handler.rfile = io.BytesIO(payload)
        handler.send_response = lambda s: status_holder.append(s)
        handler.send_header = lambda name, value: headers.__setitem__(name, value)
        handler.end_headers = lambda: None
        handler.wfile = io.BytesIO()
        handler.do_POST()
        self.assertEqual(status_holder[0], 200)
        resp = json.loads(handler.wfile.getvalue())
        self.assertEqual(resp["genomes"], ["LP-01", "LP-02"])

        # Invalid genome count (< 2)
        status_holder.clear()
        bad_payload = json.dumps({"genome_ids": ["LP-01"]}).encode()
        handler.headers = {"Content-Length": str(len(bad_payload))}
        handler.rfile = io.BytesIO(bad_payload)
        handler.wfile = io.BytesIO()
        handler.do_POST()
        self.assertEqual(status_holder[0], 400)

        # Invalid path
        status_holder.clear()
        handler.path = "/api/other"
        handler.do_POST()
        self.assertEqual(status_holder[0], 404)

    def test_real_refseq_panaroo_dataset_loading(self):
        from app import load_data, REAL_DATA
        if REAL_DATA.is_file():
            data = derive(load_data("real-lp-refseq-panaroo"))
            self.assertEqual(data["dataset"]["id"], "real-lp-refseq-panaroo")
            self.assertEqual(len(data["strains"]), 20)
            self.assertEqual(len(data["families"]), 5488)
            self.assertEqual(sum(data["counts"].values()), 5488)


if __name__ == "__main__":
    unittest.main()
