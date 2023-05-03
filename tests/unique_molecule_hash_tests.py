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
        h1 = unique_molecule_hash.get_unique_hash(mx1)
        h2 = unique_molecule_hash.get_unique_hash(mx2)
        self.assertEqual(h1, h2, "Hashes for same mixtures with different order did not match.")

    def test_absolute_stereo(self):
        """
        Test that no error happens in case a mixture contains 2 molecules with absolute stereo chemistry
        StereoGroup for absolute stereochemistry is only one group with all atoms from both components
        Therefore code needs to split them up and create a group for each component
        """
        molv3000 = """ACS Document 1996
  ChemDraw04282311222D

  0  0  0     0  0              0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 18 16 0 0 1
M  V30 BEGIN ATOM
M  V30 1 C -4.274511 0.339835 0.000000 0
M  V30 2 C -3.560040 -0.072665 0.000000 0
M  V30 3 C -2.845569 0.339835 0.000000 0
M  V30 4 C -2.131098 -0.072665 0.000000 0
M  V30 5 C -1.416627 0.339835 0.000000 0
M  V30 6 C -0.702156 -0.072665 0.000000 0
M  V30 7 O -2.845569 1.164835 0.000000 0
M  V30 8 C -4.988982 -0.072665 0.000000 0
M  V30 9 N -2.131098 -0.897665 0.000000 0
M  V30 10 C 1.416272 0.072666 0.000000 0
M  V30 11 C 2.130700 -0.339835 0.000000 0
M  V30 12 C 2.845127 0.072666 0.000000 0
M  V30 13 C 3.560128 -0.339835 0.000000 0
M  V30 14 C 4.274555 0.072666 0.000000 0
M  V30 15 C 4.988982 -0.339835 0.000000 0
M  V30 16 O 2.845127 0.897666 0.000000 0
M  V30 17 C 0.701845 -0.339835 0.000000 0
M  V30 18 N 3.560128 -1.164835 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 5 6
M  V30 6 1 3 7 CFG=1
M  V30 7 1 1 8
M  V30 8 1 4 9 CFG=1
M  V30 9 2 10 11
M  V30 10 1 11 12
M  V30 11 1 12 13
M  V30 12 1 13 14
M  V30 13 1 14 15
M  V30 14 1 12 16 CFG=3
M  V30 15 1 10 17
M  V30 16 1 13 18 CFG=1
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEABS ATOMS=(4 3 4 12 13)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
"""
        m = Chem.MolFromMolBlock(molv3000)
        h = unique_molecule_hash.get_unique_hash(m)
        self.assertIsNotNone(h)

    def test_sgroups(self):
        """
        Atoms can be annotated with text like a "*" in this case. this should not impact the hash and hash should be
        the same with or without this annotation
        """

        molv3000 = """ACS Document 1996
  ChemDraw11172216522D

  0  0  0     0  0              0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 13 13 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.681708 -1.477962 0.000000 0
M  V30 2 C -0.856798 -1.473379 0.000000 0
M  V30 3 C -0.605888 -0.687424 0.000000 0
M  V30 4 O -1.275363 -0.206227 0.000000 0
M  V30 5 C -1.940638 -0.694299 0.000000 0
M  V30 6 O -2.726593 -0.444152 0.000000 0
M  V30 7 C -0.396606 0.110753 0.000000 0
M  V30 8 C 0.192288 -0.896706 0.000000 0
M  V30 9 C 0.772780 -0.310486 0.000000 0
M  V30 10 C 1.570574 -0.519769 0.000000 0
M  V30 11 C 2.151065 0.066451 0.000000 0
M  V30 12 C 2.146483 0.891361 0.000000 0
M  V30 13 C 2.726593 1.477962 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 5 1
M  V30 6 2 5 6
M  V30 7 1 3 7
M  V30 8 1 3 8
M  V30 9 1 8 9
M  V30 10 1 9 10
M  V30 11 2 10 11
M  V30 12 1 11 12
M  V30 13 1 12 13
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 1 FIELDNAME=Text FIELDDISP="   -0.4799   -0.7485    DA    ALL  1 -
M  V30       5" FIELDDATA="*"
M  V30 END SGROUP
M  V30 END CTAB
M  END
"""
        m = Chem.MolFromMolBlock(molv3000)
        h = unique_molecule_hash.get_unique_hash(m)
        self.assertIsNotNone(h)
        m2 = Chem.MolFromSmiles("O=C1OC(C)(CC/C=C\CC)CC1")
        h2 = unique_molecule_hash.get_unique_hash(m2)
        self.assertEqual(h, h2, msg="Error: Fielddata Text changes hash.")