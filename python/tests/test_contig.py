import json,unittest
from pathlib import Path
from mixcrpy.contig import assemble_points,point_evidence,assemble_contig
from mixcrpy.reference import load_reference

class ContigTests(unittest.TestCase):
    def test_actual_original_contig_command(self):
        reference=load_reference()
        cases=json.loads(Path(__file__).with_name('contig_workflow_oracle.json').read_text())['cases']
        for case in cases:
            with self.subTest(dataset=case['dataset'],clone=case['clone']['clone_id']):
                result=assemble_contig(case['records'],case['clone'],reference)
                self.assertEqual(len(result),len(case['expected']))
                for actual,expected in zip(result,case['expected']):
                    self.assertEqual(actual['count'],expected['count'])
                    self.assertEqual([{k:t[k] for k in ('sequence','quality')} for t in actual['targets']],expected['targets'])
                    for kind in 'vdjc':
                        hits=[dict(gene=h['gene'],score=h['score'],alignments=[{k:a[k] for k in ('score','ref_start','ref_end','query_start','query_end','operations')} if a else None for a in h['alignments']]) for h in actual[kind+'_hits']]
                        self.assertEqual(hits,expected[kind+'_hits'])

    def test_original_alignment_to_contig(self):
        reference=load_reference()
        cases=json.loads(Path(__file__).with_name('point_evidence_oracle.json').read_text())['cases']
        for case in cases:
            with self.subTest(case=case['id']):
                clone=dict(cdr3_nt=case['cdr3'],count=len(case['inputs']),v_hits=[{'gene':case['v_gene']}],j_hits=[{'gene':case['j_gene']}])
                rows,position=point_evidence(case['inputs'],clone,reference)
                expected=[{int(k):tuple(v) for k,v in row.items()} for row in case['point_reads']]
                self.assertEqual(position,case['feature_position'])
                self.assertEqual(rows,expected)
                self.assertEqual(assemble_points(rows,clone['cdr3_nt'],position,clone['count']),case['expected'])

    def test_original_variant_clustering_contig_cleanup(self):
        cases=json.loads(Path(__file__).with_name('contig_oracle.json').read_text())['cases']
        for case in cases:
            with self.subTest(case=case['id']):
                rows=[{int(k):v for k,v in row.items()} for row in case['point_reads']]
                output=assemble_points(rows,case['cdr3'],case['feature_position'],case['count'],case['options'])
                self.assertEqual(len(output),len(case['expected']))
                for actual,expected in zip(output,case['expected']):
                    self.assertAlmostEqual(actual.pop('count'),expected['count'],places=10)
                    self.assertEqual(actual,{k:v for k,v in expected.items() if k!='count'})

if __name__=='__main__':unittest.main()
