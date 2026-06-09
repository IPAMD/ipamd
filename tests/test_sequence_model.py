"""Unit tests for ipamd.public.models.sequence.Sequence and subclasses."""
# pylint: disable=missing-function-docstring
# pylint: disable=missing-class-docstring
import unittest

from ipamd.public.models.sequence import (
    DNASequence,
    ProteinSequence,
    RNASequence,
    Sequence,
)


class TestSequenceModel(unittest.TestCase):
    def test_init_accepts_valid_string(self):
        seq = DNASequence("gene1", "ACGT")
        self.assertEqual(seq.name, "gene1")
        self.assertEqual(seq.sequence, "ACGT")

    def test_init_accepts_valid_list(self):
        seq = ProteinSequence("p1", ["A", "C", "G"])
        self.assertEqual(list(seq.sequence), ["A", "C", "G"])

    def test_init_rejects_non_string_non_list(self):
        with self.assertRaises(TypeError):
            DNASequence("x", 123)

    def test_init_rejects_invalid_residue(self):
        with self.assertRaises(ValueError):
            DNASequence("x", "ACGZ")

    def test_str_representation(self):
        seq = DNASequence("n1", "AC")
        self.assertEqual(str(seq), "n1:AC")

    def test_len(self):
        seq = RNASequence("r", "ACGU")
        self.assertEqual(len(seq), 4)

    def test_getitem_index_and_slice(self):
        seq = DNASequence("x", "ACGT")
        self.assertEqual(seq[0], "A")
        self.assertEqual(seq[1:3], "CG")

    def test_iter_yields_residues(self):
        seq = ProteinSequence("p", "AC")
        self.assertEqual(list(iter(seq)), ["A", "C"])

    def test_property_sequence(self):
        seq = DNASequence("n", "A")
        self.assertEqual(seq.sequence, "A")

    def test_property_name(self):
        seq = DNASequence("myname", "A")
        self.assertEqual(seq.name, "myname")

    def test_setitem_by_int_index(self):
        seq = DNASequence("x", "ACGT")
        seq[0] = "T"
        self.assertEqual(seq.sequence, "TCGT")

    def test_setitem_by_slice(self):
        seq = DNASequence("x", "ACGT")
        seq[1:3] = "TG"
        self.assertEqual(seq.sequence, "ATGT")

    def test_setitem_by_str_replace(self):
        seq = DNASequence("x", "ACAC")
        seq["AC"] = "TG"
        self.assertEqual(seq.sequence, "TGTG")

    def test_setitem_rejects_unsupported_key_type(self):
        seq = DNASequence("x", "A")
        with self.assertRaisesRegex(TypeError, "Unsupported index type"):
            seq[1.0] = "A"

    def test_setitem_rejects_invalid_value(self):
        seq = DNASequence("x", "A")
        with self.assertRaises(ValueError):
            seq[0] = "Z"

    def test_add_two_sequences_same_class(self):
        a = DNASequence("n", "AC")
        b = DNASequence("ignored", "GT")
        c = a + b
        self.assertIsInstance(c, DNASequence)
        self.assertEqual(c.name, "n")
        self.assertEqual(c.sequence, "ACGT")

    def test_add_two_sequences_different_subclass_raises(self):
        d = DNASequence("d", "A")
        r = RNASequence("r", "A")
        with self.assertRaises(TypeError):
            _ = d + r

    def test_add_string_suffix(self):
        seq = DNASequence("n", "AC")
        out = seq + "GT"
        self.assertEqual(out.sequence, "ACGT")
        self.assertEqual(out.name, "n")

    def test_add_rejects_unsupported_type(self):
        seq = DNASequence("n", "A")
        with self.assertRaises(TypeError):
            _ = seq + []

    def test_base_sequence_allows_empty_only(self):
        empty = Sequence("e", "")
        self.assertEqual(empty.sequence, "")

if __name__ == "__main__":
    unittest.main()
