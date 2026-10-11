import json
from pathlib import Path
import unittest
from mixcrpy.kalignment import extend_hit, shift_homopolymer_indels


class AlignmentExtensionTests(unittest.TestCase):
    def test_original_java_oracle(self):
        folder=Path(__file__).parent
        cases=json.loads((folder/'kmapper_oracle.json').read_text())['cases']
        answers=json.loads((folder/'kalignment_oracle.json').read_text())['cases']
        self.assertEqual(len(cases),len(answers))
        for case,expected in zip(cases,answers):
            self.assertEqual(len(case['expected']),len(expected))
            for mapper,target in zip(case['expected'],expected):
                with self.subTest(case=case['id'],gene=target['gene']):
                    target=dict(target)
                    normalized=target.pop('postprocessed_operations')
                    hit={**mapper,'gene':case['genes'][mapper['id']],'seed_positions':case['seed_positions']}
                    actual=extend_hit(case['sequence'],hit,case['start'],case['end'],case['kind'])
                    actual.pop('boundary_map',None)
                    self.assertEqual(target,actual)
                    shifted=shift_homopolymer_indels(actual,hit['gene']['sequence'],case['sequence'])
                    self.assertEqual(normalized,shifted['operations'])


if __name__=='__main__':unittest.main()
