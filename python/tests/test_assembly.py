"""Assembler behavior tests independent of any reference FASTQ fixture."""
import unittest
from mixcrpy.assembly import assemble, mapping_threshold, mapping_variants_count, mutation_positions, clonal_compatible


def alignment(sequence="TGTGCCAGCAGCTTT", quality=None, read_id=0, **extra):
    result = dict(read_id=read_id, name=str(read_id), cdr3_nt=sequence,
                  cdr3_aa="CASSF", cdr3_quality=quality or "I"*len(sequence),
                  v_hits=[{"gene": "TRBV1*00", "score": 200}],
                  j_hits=[{"gene": "TRBJ1*00", "score": 100}],
                  c_hits=[], d_hits=[], n_regions=[[6, 9]])
    result.update(extra)
    return result


class AssemblyTests(unittest.TestCase):
    def test_exact_grouping_score_and_max_quality(self):
        a = alignment(quality="H"*15)
        b = alignment(read_id=1, v_hits=[{"gene": "TRBV1*00", "score": 150}])
        result = assemble([a, b])
        clone = result["clones"][0]
        self.assertEqual((clone["count"], clone["fraction"], clone["quality"]), (2, 1, "I"*15))
        self.assertEqual(clone["v_hits"][0]["score"], 175)
        self.assertEqual(result["report"]["assigned_reads"], 2)

    def test_low_quality_only_mutates_low_quality_positions(self):
        core = alignment()
        near = alignment(sequence="TGTGCCAGCAGATTT", quality="I"*11+"!"+"I"*3, read_id=1)
        bad = alignment(sequence="TGTGCCAGCAGATTT", quality="!"+"I"*14, read_id=2)
        result = assemble([core, near, bad])
        self.assertEqual(result["clones"][0]["count"], 2)
        self.assertEqual(result["clones"][0]["mapped_count"], 1)
        self.assertEqual(result["clones"][0]["quality"], "I"*15)
        self.assertEqual(result["read_assignments"][2]["reason"], "no_candidate_for_deferred_alignment")

    def test_q20_is_bad_and_fraction_filter_is_strict(self):
        result = assemble([alignment(quality="5"*15)])
        self.assertEqual(result["report"]["too_many_low_quality_points"], 1)
        self.assertEqual(result["clones"], [])

    def test_clustering_read_counts_default_and_option(self):
        reads = [alignment(read_id=i) for i in range(100)]
        reads.append(alignment(sequence="TGTGCCAGCAGATTT", read_id=100))
        result = assemble(reads)
        self.assertEqual(len(result["clones"]), 1)
        self.assertEqual(result["clones"][0]["count"], 100)
        self.assertEqual(result["report"]["assigned_reads"], 101)
        self.assertEqual(result["report"]["clustered_reads"], 1)
        self.assertTrue(result["read_assignments"][-1]["clustered"])
        added = assemble(reads, {"addReadsCountOnClustering": True})
        self.assertEqual(added["clones"][0]["count"], 101)

    def test_relative_concentration_retains_supported_minor(self):
        reads = [alignment(read_id=i) for i in range(20)]
        reads.append(alignment(sequence="TGTGCCAGCAGATTT", read_id=20))
        self.assertEqual(len(assemble(reads)["clones"]), 2)

    def test_n_regions_limit_counts_affected_regions(self):
        reads = [alignment(read_id=i, n_regions=[[5, 7], [8, 11]]) for i in range(5000)]
        reads.append(alignment(sequence="TGTGCCCGCAGATTT", read_id=5000))
        self.assertEqual(len(assemble(reads)["clones"]), 2)

    def test_edit_positions_and_indels(self):
        self.assertEqual(mutation_positions("ACTG", "ATTG"), (1,))
        self.assertEqual(len(mutation_positions("ACTG", "ACTTG")), 1)
        self.assertIsNone(mutation_positions("ACTG", "AAAAAAA"))

    def test_mapping_threshold_boundaries(self):
        self.assertEqual(mapping_variants_count(5, 2), 106)
        self.assertEqual([mapping_threshold(i) for i in [0, 1, 2, 5, 6, 35, 36]], [0, 1, 2, 2, 1, 1, 0])

    def test_short_cdr3_and_missing_quality(self):
        self.assertEqual(assemble([alignment(sequence="TGT")])["report"]["dropped_reads"], 1)
        with self.assertRaises(ValueError):
            assemble([alignment(cdr3_quality=None)])

    def test_invalid_options_fail_even_for_empty_input(self):
        for options in ({"bogus": 1}, {"relativeMinScores": {"x": .9}},
                        {"cloneClusteringParameters": {"unknown": True}},
                        {"cloneClusteringParameters": {"clusteringFilter": {"unknown": True}}}):
            with self.assertRaises(ValueError):
                assemble([], options)
        for options in ({"cloneClusteringParameters": {"searchParameters": "oneMismatch"}},
                        {"cloneClusteringParameters": {"clusteringFilter": {"type": "absolute"}}}):
            with self.assertRaises(NotImplementedError):
                assemble([], options)

    def test_valid_partial_options_and_disabled_clustering(self):
        self.assertEqual(assemble([alignment()], {"relativeMinScores": {"v": .8}})["clones"][0]["count"], 1)
        self.assertEqual(assemble([alignment()], {"cloneClusteringParameters": None})["clones"][0]["count"], 1)

    def test_preclustering_merges_counts_preserves_major_quality_and_scores(self):
        reads = [alignment(read_id=i, quality="H"*15) for i in range(2)]
        reads.append(alignment(read_id=2, v_hits=[{"gene":"TRBV2*00","score":205},
                                                 {"gene":"TRBV1*00","score":200}]))
        result=assemble(reads,{"separateByV":True})
        clone=result["clones"][0]
        self.assertEqual(clone["count"],3)
        self.assertEqual(clone["quality"],"H"*15)
        self.assertEqual(clone["v_hits"],[{"gene":"TRBV1*00","score":200}])
        self.assertTrue(result["read_assignments"][2]["preclustered"])
        self.assertEqual(result["report"]["reads_preclustered"],1)

    def test_multiple_features_preserve_target_boundaries(self):
        a=alignment(read_id=0,features={"A":{"sequence":"AAAAAA","quality":"I"*6},
                                       "B":{"sequence":"CCCCCC","quality":"H"*6}})
        b=alignment(read_id=1,features={"A":{"sequence":"AAAAAAC","quality":"I"*7},
                                       "B":{"sequence":"CCCCC","quality":"I"*5}})
        result=assemble([a,b],{"assemblingFeatures":["A","B"],"cloneClusteringParameters":None})
        self.assertEqual(len(result["clones"]),2)
        self.assertEqual(result["clones"][0]["clonal_sequence"],"AAAAAACCCCCC")
        self.assertEqual(result["clones"][0]["cdr3_nt"],"")
        self.assertNotEqual([len(t["sequence"]) for t in result["clones"][0]["targets"]],
                            [len(t["sequence"]) for t in result["clones"][1]["targets"]])

    def test_gap_boundary_compatibility(self):
        self.assertFalse(clonal_compatible(["AAAA","CCCC"],["AAAAA","CCC"]))
        self.assertFalse(clonal_compatible(["AAAA","CCCC"],["AAAAA","CCC"],[]))
        self.assertTrue(clonal_compatible(["AAAA","AAAA"],["AAAAA","AAAA"],[(4,"I","A")]))

    def test_partial_read_ids_expand_assignments_and_weight_counts(self):
        a=alignment(read_id=10,read_ids=[10,12],number_of_reads=2)
        b=alignment(read_id=15)
        result=assemble([a,b])
        self.assertEqual([a["read_id"] for a in result["read_assignments"]],[10,12,15])
        self.assertEqual(result["clones"][0]["count"],3)
        self.assertEqual(result["clones"][0]["v_hits"][0]["score"],133.3)
        self.assertEqual(result["report"]["assigned_reads"],3)
        self.assertEqual(result["report"]["total_alignments"],2)
        with self.assertRaises(ValueError):
            assemble([a,alignment(read_id=12)])


if __name__ == "__main__":
    unittest.main()
