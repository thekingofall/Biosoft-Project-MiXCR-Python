import unittest
from types import SimpleNamespace
from unittest.mock import patch
from mixcrpy.features import extract_feature, relative_range, target_alignment


class FeatureRangeTests(unittest.TestCase):
    def setUp(self):
        self.reference = SimpleNamespace(features={'CDR3': {
            'begin': 'CDR3Begin', 'end': 'CDR3End', 'contains_junction': True,
            'contains_d': True, 'contains_cdr3': True}})
        self.record = {'targets': [{'sequence': 'AACCGGTTAA', 'quality': 'I'*10},
                                   {'sequence': 'TTGGCCAATT', 'quality': '!'*10}],
                       'd_hits': [{'gene': 'example'}]}

    def extract(self, ranges):
        with patch('mixcrpy.features.feature_range',
                   side_effect=lambda record, ref, name, index=0: ranges.get((index, name))):
            return extract_feature(self.record, self.reference, 'CDR3')

    def test_best_quality_sequence_and_last_available_relative_range(self):
        value = self.extract({(0, 'CDR3'): (1, 9), (1, 'CDR3'): (0, 10),
                              (0, 'VJJunction'): (3, 4), (1, 'VJJunction'): (3, 8)})
        self.assertEqual(value['sequence'], 'ACCGGTTA')
        self.assertEqual(value['target_index'], 0)
        self.assertEqual(value['n_regions'], [[3, 8]])

    def test_zero_length_d_junction_blocks_vj_fallback(self):
        value = self.extract({(0, 'CDR3'): (1, 9), (0, 'VDJunction'): (4, 4),
                              (0, 'VJJunction'): (3, 8)})
        self.assertEqual(value['n_regions'], [[3, 3]])

    def test_partial_junction_is_not_clipped_and_d_hit_does_not_block_fallback(self):
        value = self.extract({(0, 'CDR3'): (1, 9), (0, 'VDJunction'): (0, 4),
                              (0, 'DJJunction'): (7, 10), (0, 'VJJunction'): (4, 4)})
        self.assertEqual(value['n_regions'], [[3, 3]])

    def test_empty_explicit_alignment_list(self):
        self.assertIsNone(target_alignment({'alignments': []}, 0))


if __name__ == '__main__':
    unittest.main()
