"""Numeric and coordinate checks independent of MiXCR fixture outputs."""
import random
import unittest
from mixcrpy.alignment import linear_align, translate, reverse_complement, _batch_scores


def scalar_score(reference, query, left, right):
    rows = [[0]*(len(query)+1) for _ in range(len(reference)+1)]
    for i,a in enumerate(reference,1):
        for j,b in enumerate(query,1):
            candidates=[rows[i-1][j-1]+(5 if a==b else -9),
                        rows[i-1][j]-12, rows[i][j-1]-12]
            if left: candidates.append(0)
            rows[i][j]=max(candidates)
    return max(max(row) for row in rows) if right else max(max(rows[-1]),max(row[-1] for row in rows))


class AlignmentTests(unittest.TestCase):
    def test_vectorized_dp_against_scalar(self):
        rng=random.Random(813713)
        for _ in range(150):
            reference=''.join(rng.choices('ACGT',k=rng.randrange(1,40)))
            query=''.join(rng.choices('ACGT',k=rng.randrange(1,40)))
            for left in (False,True):
                for right in (False,True):
                    alignment=linear_align(reference,query,floating_left=left,floating_right=right)
                    score=alignment['score'] if alignment else 0
                    self.assertEqual(score,scalar_score(reference,query,left,right))
                    if alignment:
                        ops=alignment['operations']
                        self.assertEqual(score,5*ops.count('M')-9*ops.count('X')-12*(ops.count('I')+ops.count('D')))
                        self.assertEqual(alignment['ref_end']-alignment['ref_start'],sum(op!='I' for op in ops))
                        self.assertEqual(alignment['query_end']-alignment['query_start'],sum(op!='D' for op in ops))

    def test_reverse_complement_ambiguity(self):
        sequence='ACGTRYKMSWBDHVN'
        self.assertEqual(reverse_complement(reverse_complement(sequence)),sequence)

    def test_batch_scores_against_scalar(self):
        rng=random.Random(88863)
        genes=[''.join(rng.choices('ACGT',k=rng.randrange(1,50))) for _ in range(40)]
        query=''.join(rng.choices('ACGT',k=39))
        for left in (False,True):
            for right in (False,True):
                self.assertEqual(_batch_scores(genes,query,left,right),
                                 [scalar_score(g,query,left,right) for g in genes])

    def test_center_frame_translation(self):
        self.assertEqual(translate('TGTTTT'),'CF')
        self.assertEqual(translate('TGTACTTT'),'C_F')
        self.assertEqual(translate('TGTAATTT'),'C_F')
        self.assertEqual(translate(''),'')


if __name__=='__main__': unittest.main()
