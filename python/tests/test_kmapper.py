import json
from pathlib import Path
import unittest
from mixcrpy.kmapper import map_candidates
from mixcrpy.randomutil import Well19937c


class ReplayRandom:
    def __init__(self,draws):
        self.draws=iter(draws)
    def nextInt(self,bound):
        assert bound==7
        return next(self.draws)


class MapperTests(unittest.TestCase):
    def test_original_java_oracle(self):
        data=json.loads(Path(__file__).with_name('kmapper_oracle.json').read_text())
        for case in data['cases']:
            with self.subTest(case=case['id']):
                rng=ReplayRandom(case['draws'])
                hits=map_candidates(case['sequence'],case['genes'],rng,case['start'],case['end'],case['kind'])
                with self.assertRaises(StopIteration):
                    next(rng.draws)
                values=[]
                for hit in hits:
                    self.assertEqual(case['seed_positions'],hit['seed_positions'])
                    values.append({k:v for k,v in hit.items() if k not in {'gene','seed_positions'}})
                self.assertEqual(case['expected'],values)

    def test_python_well_random_stream(self):
        cases=json.loads(Path(__file__).with_name('kmapper_oracle.json').read_text())['cases']
        for case in cases:
            with self.subTest(case=case['id']):
                rng=Well19937c(case['seed'])
                self.assertEqual(case['draws'],[rng.next_int(7) for _ in case['draws']])

    def test_empty_filter_still_consumes_stream(self):
        case=json.loads(Path(__file__).with_name('kmapper_oracle.json').read_text())['cases'][1]
        rng=ReplayRandom(case['draws'])
        self.assertEqual([],map_candidates(case['sequence'],case['genes'],rng,case['start'],case['end'],case['kind'],allowed_names=set()))
        with self.assertRaises(StopIteration):next(rng.draws)

    def test_short_window_consumes_no_random(self):
        self.assertEqual([],map_candidates('ACGTA',[{'name':'g','sequence':'ACGTACGT'}],ReplayRandom([])))


if __name__=='__main__': unittest.main()
