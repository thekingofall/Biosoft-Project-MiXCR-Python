"""Original-JAR paired alignment plus multi-target partial assembly oracle."""
import json,unittest
from pathlib import Path
from mixcrpy.reference import load_reference
from mixcrpy.partial import assemble_partial
from mixcrpy.paired import align_pair


def scientific(record):
    targets=tuple((t['sequence'],t['quality']) for t in record['targets'])
    hits=[]
    for kind in 'vdjc':
        hits.append(tuple((h['gene'],h['score'],tuple(None if a is None else tuple(a[k] for k in ('score','ref_start','ref_end','query_start','query_end','operations')) for a in h['alignments'])) for h in record[kind+'_hits']))
    return tuple(sorted(record['read_ids'])),targets,record['cdr3_nt'],record['cdr3_quality'],tuple(hits)


class PartialMultiTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cases=json.loads(Path(__file__).with_name('partial_multi_oracle.json').read_text())['cases']
        cls.reference=load_reference()

    def test_original_multitarget_assembly(self):
        for case in self.cases:
            with self.subTest(case=case['id']):
                result=assemble_partial(case['inputs'],self.reference,{'chains':case['chains']})
                self.assertEqual(sorted(map(scientific,case['outputs'])),sorted(map(scientific,result['alignments'])))
                for key,value in case['report'].items():self.assertEqual(value,result['report'].get(key,0),key)

    def test_paired_fastq_through_partial(self):
        for case in self.cases:
            with self.subTest(case=case['id']):
                records=[]
                for i,read in enumerate(case['reads']):
                    result=align_pair(i,f'r{i}/1',read['sequence1'],read['quality1'],f'r{i}/2',read['sequence2'],read['quality2'],self.reference,case['chains'],{'no_merge':True,'allowPartialAlignments':True})
                    if result['status']=='aligned':
                        result['read_ids']=[i];result['alignment_index']=len(records);records.append(result)
                self.assertEqual(list(map(scientific,case['inputs'])),list(map(scientific,records)))
                result=assemble_partial(records,self.reference,{'chains':case['chains']})
                self.assertEqual(sorted(map(scientific,case['outputs'])),sorted(map(scientific,result['alignments'])))

if __name__=='__main__':unittest.main()
