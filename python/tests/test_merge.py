import json
from pathlib import Path
import unittest
from mixcrpy.merge import merge_pair, reverse_complement


class MergeTests(unittest.TestCase):
    def test_perfect_overlap(self):
        sequence = "CGCACAGTGTTGTCAAAGAAAACGCGTACGACATTGAGAAGACCGGCCGTTCTCCTTTGACATGATTGGATCGGTTGCTGCCGGCCCAGAATCCTAGCAG"
        first, second = sequence[:70], sequence[45:]
        result = merge_pair(first, "I" * len(first), reverse_complement(second), "H" * len(second))
        self.assertEqual(sequence, result["sequence"])
        self.assertEqual(25, result["overlap"])
        self.assertEqual(45, result["offset"])

    def test_quality_disagreement(self):
        sequence = "ACGCTGACTGACTGCGTAGCTAGTGCATGCTGAC"
        altered = sequence[:10] + "T" + sequence[11:]
        result = merge_pair(sequence, "I" * len(sequence), reverse_complement(altered), "5" * len(sequence))
        self.assertEqual(sequence, result["sequence"])
        self.assertEqual("5", result["quality"][10])
        self.assertEqual(1, result["errors"])

    def test_too_short_and_malformed(self):
        self.assertIsNone(merge_pair("A" * 16, "I" * 16, "T" * 16, "I" * 16))
        with self.assertRaises(ValueError):
            merge_pair("AA", "I", "TT", "II")
        with self.assertRaises(ValueError):
            merge_pair("N" * 20, "I" * 20, "N" * 20, "I" * 20)

    def test_original_java_oracle(self):
        path = Path(__file__).with_name("merge_oracle.json")
        data = json.loads(path.read_text())
        for case in data["cases"]:
            with self.subTest(case=case["id"]):
                self.assertEqual(case["expected"], merge_pair(*case["input"]))


if __name__ == "__main__":
    unittest.main()
