from unique_molecule_hash import *
import unittest
from rdkit import Chem
from pathlib import Path
import logging

logger = logging.getLogger('unique_molecule_hash')
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.DEBUG)


class UniqueMoleculeHashTest(unittest.TestCase):
    """
    Test if conversions match expected outcome
    """

    def test_mixture_order(self):
        """
        Test that component order does not matter
        """
        mx1 = Chem.MolFromSmiles("CC.CCC.CC=O")
        mx2 = Chem.MolFromSmiles("CCC.CC=O.CC")
        h1 = unique_molecule_hash.get_hash(mx1)
        h2 = unique_molecule_hash.get_hash(mx2)
        self.assertEqual(h1, h2, "Hashes for same mixtures with different order did not match.")

    def test_absolute_stereo(self):
        """
        Test that no error happens in case a mixture contains 2 molecules with absolute stereo chemistry
        StereoGroup for absolute stereochemistry is only one group with all atoms from both components
        Therefore code needs to split them up and create a group for each component
        """
        m = Chem.MolFromMolFile(r"files/mixture_of_isomers.mol")
        h = unique_molecule_hash.get_hash(m)
        self.assertIsNotNone(h)

    def test_sgroups(self):
        """
        Atoms can be annotated with text like a "*" in this case. this should not impact the hash and hash should be
        the same with or without this annotation
        """

        m = Chem.MolFromMolFile(r"files/field_data_txt.mol")
        h = unique_molecule_hash.get_hash(m)
        self.assertIsNotNone(h)
        m2 = Chem.MolFromSmiles("O=C1OC(C)(CC/C=C\CC)CC1")
        h2 = unique_molecule_hash.get_hash(m2)
        self.assertEqual(h, h2, msg="Error: Fielddata Text changes hash.")

    def test_brackets(self):

        pm1 = Chem.MolFromMolFile(r"files/polymer1.mol")
        pm2 = Chem.MolFromMolFile(r"files/polymer2.mol")

        h1 = unique_molecule_hash.get_hash(pm1)
        h2 = unique_molecule_hash.get_hash(pm2)
        self.assertIsNotNone(h1)
        self.assertIsNotNone(h2)
        self.assertIsNot(h1,h2,"Polymers have same hash even-though different repeat pattern specified.")

    def test_tautomer_insensitivity(self):

        smi = "CC/C=C1/C2C(=CC=CC=2)C(=O)O1"
        m = Chem.MolFromSmiles(smi)

        # Generate random kekulized SMILES and check if hash is always the same
        hashes = []
        for i in range(100):
            hashes.append(unique_molecule_hash.get_standard_hash(
                Chem.MolFromSmiles(Chem.MolToSmiles(m, doRandom=True, canonical=False, kekuleSmiles=True))))

        self.assertTrue(len(set(hashes)) == 1,
                        "Found random tautomer with different hash. Hash is not tautomer insensitive")

    def test_enhanced_stereo_input(self):

        or1 = Chem.MolFromSmiles("C[C@@H](O)CC |o1:1|")
        h_or1 = unique_molecule_hash.get_standard_hash(or1)
        or2 = Chem.MolFromSmiles("C[C@H](O)CC |o1:1|")
        h_or2 = unique_molecule_hash.get_standard_hash(or2)
        self.assertEqual(h_or1, h_or2, "Enhanced Stereo 'OR' with different bond direction do not match.")

        and1 = Chem.MolFromSmiles("C[C@@H](O)CC |&1:1|")
        h_and1 = unique_molecule_hash.get_standard_hash(and1)
        and2 = Chem.MolFromSmiles("C[C@H](O)CC |&1:1|")
        h_and2 = unique_molecule_hash.get_standard_hash(and2)
        self.assertEqual(h_and1, h_and2, "Enhanced Stereo 'AND' with different bond direction do not match.")

        abs = Chem.MolFromSmiles("C[C@@H](O)CC |a:1|")
        h_abs = unique_molecule_hash.get_standard_hash(abs)
        self.assertNotEqual(h_abs, h_or1, "Molecule with absolute stereo has same hash as with 'OR'.")
        self.assertNotEqual(h_abs, h_and1, "Molecule with absolute stereo has same hash as with 'AND'.")

    def test_attachment_point(self):

        m1 = Chem.MolFromMolFile(r"files/attachment_point1.mol")
        m2 = Chem.MolFromMolFile(r"files/attachment_point2.mol")

        h1 = unique_molecule_hash.get_hash(m1)
        h2 = unique_molecule_hash.get_hash(m2)
        self.assertIsNotNone(h1)
        self.assertIsNotNone(h2)
        self.assertEqual(h1, h2, "Attachment point molecules have a different hash. ")
