import json
from pathlib import Path
import unittest
from mixcrpy.partial import assemble_partial
from mixcrpy.reference import load_reference


class PartialTests(unittest.TestCase):
    def test_full_original_partial_assembler(self):
        reference=load_reference()
        cases=json.loads(Path(__file__).with_name('partial_oracle.json').read_text())['cases']
        for case in cases:
            for record in case['inputs']:
                for kind in 'vdjc':
                    for hit in record[kind+'_hits']:
                        hit['boundary_map']={int(k):v for k,v in hit['boundary_map'].items()}
            with self.subTest(case=case['id']):
                result=assemble_partial(case['inputs'],reference,{'chains':case['chains']})
                for key,value in case['report'].items():
                    self.assertEqual(value,result['report'].get(key,0),key)
                def scientific(record):
                    return (tuple(record['read_ids']),record['sequence'],record['quality'],record['cdr3_nt'],record['cdr3_quality'],
                            tuple((kind,tuple((h['gene'],h['score'],h['ref_start'],h['ref_end'],h['query_start'],h['query_end']) for h in record[kind+'_hits'])) for kind in 'vdjc'))
                self.assertEqual(sorted(map(scientific,case['outputs'])),sorted(map(scientific,result['alignments'])))


if __name__=='__main__':unittest.main()
