"""Fixed original-JAR outputs across paired and partial modes."""
import json,re,unittest
from pathlib import Path
from mixcrpy.reference import load_reference
from mixcrpy.alignment import align_read
from mixcrpy.paired import align_pair

class PairedOracleTests(unittest.TestCase):
    def test_original_jar_outputs(self):
        refs={f:load_reference(v_feature=f) for f in ('VRegionWithP','VGeneWithP')}
        cases=json.loads(Path(__file__).with_name('paired_oracle.json').read_text())['cases']
        for case in cases:
            with self.subTest(dataset=case['dataset'],read_id=case['read_id']):
                reads=case['reads'];args=[]
                for read in reads:args.extend(read[k] for k in ('name','sequence','quality'))
                ref=refs[case['v_feature']]
                result=(align_pair if len(reads)==2 else align_read)(case['read_id'],*args,ref,options=case['options'])
                expected=case['expected']
                self.assertEqual(bool(expected),result['status']=='aligned')
                if expected is None:continue
                for kind in 'VDJC':
                    want=[(a,float(b)) for a,b in re.findall(r'([^,()]+)\(([\d.]+)\)',expected['all'+kind+'HitsWithScore'])]
                    actual=[(h['gene'],float(h['score'])) for h in result[kind.lower()+'_hits']]
                    self.assertEqual(want,actual)
                for src,dst in [('nSeqCDR3','cdr3_nt'),('aaSeqCDR3','cdr3_aa'),('qualCDR3','cdr3_quality')]:
                    self.assertEqual(expected[src],result[dst])
                targets=result.get('targets',[dict(sequence=result['oriented_sequence'],quality=result['oriented_quality'])])
                self.assertEqual(expected['targetSequences'],','.join(t['sequence'] for t in targets))
                self.assertEqual(expected['targetQualities'],','.join(t['quality'] for t in targets))

    def test_unknown_options_rejected(self):
        ref=load_reference()
        with self.assertRaisesRegex(ValueError,'Unsupported alignment options'):
            align_read(0,'r','ACGT','IIII',ref,options={'typo':True})
        with self.assertRaisesRegex(ValueError,'Unsupported alignment options'):
            align_pair(0,'r','ACGT','IIII','r','ACGT','IIII',ref,options={'typo':True})

if __name__=='__main__':unittest.main()
