import unittest
from mixcrpy.randomutil import Well19937c, read_seed, replace_wildcards, java_string_hash


class RandomUtilTests(unittest.TestCase):
    def test_java_bounded_integer_stream(self):
        generator = Well19937c(0)
        bounds = [1,2,3,4,7,16,127,1024,10000,1073741825,2147483647]
        expected = [0,1,0,0,2,3,4,396,2178,874745425,395315341]
        self.assertEqual([generator.nextInt(n) for n in bounds], expected)

    def test_java_long_stream_and_reseed(self):
        generator = Well19937c(0)
        expected = [-4736355256857903122, -7529656944096573677, 6601006367108760060]
        self.assertEqual([generator.nextLong() for _ in range(3)], expected)
        generator.reseed(0)
        self.assertEqual(generator.nextLong(), expected[0])

    def test_read_seed_description_and_quality(self):
        self.assertEqual(read_seed(0,"read 1 /1","AGCTAGCTAGCT"), -46329256949)
        self.assertEqual(read_seed(9,"read 1 /1","AGCTAGCTAGCT","!"*12), -46329256949)
        self.assertNotEqual(read_seed(0,"read","AGCT"), read_seed(0,"read extra","AGCT"))

    def test_iupac_replacement_precedes_hash(self):
        sequence = "ACGTNRYSWKMBDHVN"
        self.assertEqual(replace_wildcards(0,sequence,"I"*len(sequence)),
                         ("ACGTAGTGAGATACGG", "IIII!!!!!!!!!!!!"))
        self.assertEqual(read_seed(0,"read 1 /1",sequence), 26279545288)

    def test_java_string_utf16_and_bound_validation(self):
        self.assertEqual(java_string_hash("\U0001f600"), 1772899)
        for bound in (0, -1, 2147483648):
            with self.assertRaises(ValueError):
                Well19937c(0).next_int(bound)


if __name__ == "__main__":
    unittest.main()
