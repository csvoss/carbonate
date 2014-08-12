"""
Unit Tests for toMolecule.py
"""

import itertools
from random import random, randrange
import unittest

from molecularStructure import Atom, Molecule
from toMolecule import moleculify, BASIC_BONDS
from toSmiles import smilesify, to_canonical


DEBUG = True

######################
##### UNIT TESTS #####
######################

def example_molecule():
    """
    This method intended as a code example of how to create and return
    a molecule.
    This is useful to reference while writing the parser.
    return :: Molecule.
    """
    c40 = Atom("C")
    c41 = Atom("C")
    mol4 = Molecule(c40)
    mol4.addAtom(c41, c40, 2)
    c42 = Atom("C")
    c43 = Atom("C")
    mol4.addAtom(c42, c40, 1)
    mol4.addAtom(c43, c41, 1)
    cl1 = Atom("Cl")
    cl2 = Atom("Cl")
    mol4.addAtom(cl1, c40, 1)
    mol4.addAtom(cl2, c41, 1)
    c40.newCTCenter(c41, c42, cl1)
    c41.newCTCenter(c40, c43, cl2)
    mol4.addBond(c42, c43, 1)

    return mol4.withHydrogens()

example_smiles_easy = [
    r"C1CCCCC1",         ## cyclohexane
    r"N#N",              ## dinitrogen
    r"CN=C=O",           ## methyl isocyanate
    r"N1CCN(CC1)C(C(F)=C2)=CC(=C2C4=O)N(C3CC3)C=C4C(=O)O",
    r"C12C3(CO)C4C1(CC(=O)O)C5(N)C2(CC(=O)C)C3(C=CC=C)C45",
    r"C(C=C2)C1C(=O)C(C(C)(C)C)C(CN)CC21",
    r"N1CCN(CC1)C(C(F)=C2)=CC(=C2C4=O)N(C3CC3)C=C4C(=O)O",
    r"C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C))))))))))))))))))))C",
]
example_smiles_atoms = [
    r"OC(=O)CC(P)C(SBr)C(CB)C(F)(Cl)C(I)Br",
    r"CC[*]CCC[*]CC",
    r"[2H]C(Cl)(Cl)Cl",  ## deuterochloroform (hydrogen-2)
    r"[Cu+2].[O-]S(=O)(=O)[O-]" ## copper(ii) sulfate
    r"C(CC[Se]C)C(CC(O)CC)Br",
    r"[Rh-](Cl)(Cl)(Cl)(Cl)$[Rh-](Cl)(Cl)(Cl)Cl",
    r'C(CCC)C([NH7])Br',
    r'P(=O)([O-])([O-])[O-]',
    r'C(=O)([O-])[O-]',
    r'C12346F5O1Cl2Br3I4S65',
    r'[Cl-].[Na+]',
    r'C#[C-]',
]
example_smiles_cistrans = [
    r"C(/F)(\Cl)=C(/Br)\I",  ## cis/trans test
    r"F\C=C(/Br)\I",     ## more cis/trans
    r"F\C=C/Br",         ## more cis/trans
    r"F/C=C/F",          ## trans-difluoroethene
    r"F/C=C\F",          ## cis-difluoroethene
]
example_smiles_chiral = [
    r"N[C@@H](C)C(=O)O", ## L-alanine
    r"N[C@H](C(=O)O)C",  ## also L-alanine
    r"N[C@H](C)C(=O)O",  ## D-alanine (less common)
    r"C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O", ## glucose
    r"CC[C@H](O1)CC[C@@]12CCCO2", ## another pheromone
    r"CC(C)[C@@]12C[C@@H]1[C@@H](C)C(=O)C2", ## alpha-thujone
    r"CC([C@H](C)Cl)Br",
    r"C([C@H](F)Cl)",
]
example_smiles_hard = [
    r"CCC[C@@H](O)CC\C=C\C=C\C#CC#C\C=C\CO", ## enanthotoxin
    r"COC(=O)C(\C)=C\C1C(C)(C)[C@H]1C(=O)O[C@@H]2C(C)=C(C(=O)C2)CC="+\
        r"CC=C", ## pyrethrin II
    r"CC(=O)OCCC(/C)=C\C[C@H](C(C)=C)CCC=C", ## some pheromone
    r"C[C@@](C)(O1)C[C@@H](O)[C@@]1(O2)[C@@H](C)[C@@H]3CC=C4[C@]3(C"+\
        "2)C(=O)C[C@H]5[C@H]4CC[C@@H](C6)[C@]5(C)CC(N7)C6NC(C[C@@]"+\
        "89(C))C7C[C@@H]8CC[C@@H]%10[C@@H]9C[C@@H](O)[C@@]%11(C)C%"+\
        "10=C[C@H](O%12)[C@]%11(O)[C@H](C)[C@]%12(O%13)[C@H](O)C[C"+\
        "@@]%13(C)CO", ## Cephalostatin-1 without aromaticity. Note
                       ## the % before ring closure labels above 9.
]
example_smiles_hard_no_cistrans = [
    r"CCC[C@@H](O)CCC=CC=CC#CC#CC=CCO", ## enanthotoxin
    r"COC(=O)C(C)=CC1C(C)(C)[C@H]1C(=O)O[C@@H]2C(C)=C(C(=O)C2)CC=CC=C",
    r"CC(=O)OCCC(C)=CC[C@H](C(C)=C)CCC=C", ## some pheromone
    r"C[C@@](C)(O1)C[C@@H](O)[C@@]1(O2)[C@@H](C)[C@@H]3CC=C4[C@]3(C"+\
        "2)C(=O)C[C@H]5[C@H]4CC[C@@H](C6)[C@]5(C)CC(N7)C6NC(C[C@@]"+\
        "89(C))C7C[C@@H]8CC[C@@H]%10[C@@H]9C[C@@H](O)[C@@]%11(C)C%"+\
        "10=C[C@H](O%12)[C@]%11(O)[C@H](C)[C@]%12(O%13)[C@H](O)C[C"+\
        "@@]%13(C)CO", ## Cephalostatin-1 without aromaticity. Note
                       ## the % before ring closure labels above 9.
]

example_smiles_hard_no_chirality_or_cistrans = [
    r"CCC[CH](O)CCC=CC=CC#CC#CC=CCO", ## enanthotoxin
    r"COC(=O)C(C)=CC1C(C)(C)[CH]1C(=O)O[CH]2C(C)=C(C(=O)C2)CC=CC=C",
    r"CC(=O)OCCC(C)=CC[CH](C(C)=C)CCC=C", ## some pheromone
    r"C[C](C)(O1)C[CH](O)[C]1(O2)[CH](C)[CH]3CC=C4[C]3(C2)C(=O)C[CH]5[CH]4CC[CH](C6)[C]5(C)CC(N7)C6NC(C[C]89(C))C7C[CH]8CC[CH]%10[CH]9C[CH](O)[C]%11(C)C%10=C[CH](O%12)[C]%11(O)[CH](C)[C]%12(O%13)[CH](O)C[C]%13(C)CO", ## Cephalostatin-1 without aromaticity. Note the % before ring closure labels above 9.
]
# C[C](C)(O1)C[CH](O)[C]1(O2)[CH](C)[CH]3CC=C4[C]3(C2)C(=O)C[CH]5[CH]4CC[CH](C6)[C]5(C)CC(N7)C6NC(C[C]89(C))C7C[CH]8CC[CH]%2510[CH]9C[CH](O)[C]%2511(C)C%2510=C[CH](O%2512)[C]%2511(O)[CH](C)[C]%2512(O%2513)[CH](O)C[C]%2513(C)CO
# OCC1(C)CC(C2(O1)OC1C(C2C)(O)C2(C(=C1)C1CCC3C(C1CC2O)(C)CC1C(C3)NC2C(N1)CC1C(C2)(C)C2CC(=O)C34C(=CCC4C(C4(OC3)OC(CC4O)(C)C)C)C2CC1)C)O

example_smiles = example_smiles_easy + example_smiles_cistrans +\
    example_smiles_chiral + example_smiles_hard + example_smiles_atoms


class TestToMolecule(unittest.TestCase):
    "Test ALL the features!"

    def assertOne(self, smi):
        "Assert that a single smiles is OK."
        try:
            self.assertEqual(
                to_canonical(smi),
                smilesify(moleculify(smi), canonical=True)
            )
            self.assertEqual(
                to_canonical(smi),
                smilesify(moleculify(to_canonical(smi)), canonical=True)
            )
        except AssertionError, StandardError:
            if DEBUG:
                print "\n%s\n" % smi
                raise

    def assertMany(self, smileses):
        "Assert that each in a list of smiles are OK."
        for smi in smileses:
            self.assertOne(smi)

    def assertSame(self, smiles1, smiles2):
        "Assert that two smiles produce the same output."
        self.assertEqual(
            smilesify(moleculify(smiles1), canonical=True),
            smilesify(moleculify(smiles2), canonical=True)
        )
        self.assertEqual(
            to_canonical(smiles1),
            smilesify(moleculify(smiles1), canonical=True)
        )

    def assertSameMany(self, smileses):
        "Assert that each in a list of smiles are identical."
        for combo in itertools.combinations(smileses, 2):
            self.assertSame(combo[0], combo[1])

    def test_carbons(self):
        "Test some basic strings of carbon."
        smileses = ["C"*num for num in range(1, 21)]
        self.assertMany(smileses)

    def test_branches(self):
        "Test basic tree-branching carbon strings. Randomly-generated."
        termination_probability = 0.90
        max_run_length = 5
        def random_branching_carbons():
            "Create a randomly-branching carbon string."
            def end_or_branch():
                "Terminate or branch further."
                if random() < termination_probability:
                    return "C" * (1 + randrange(max_run_length))
                else:
                    return random_branching_carbons()
            a = end_or_branch()
            b = end_or_branch()
            c = end_or_branch()
            d = end_or_branch()
            if random() < 0.5:
                return a + "(" + b + ")(" + c + ")" + d
            else:
                return a + "(" + b + ")" + c

        for _ in range(1, 20):
            try:
                smi = random_branching_carbons()
            except RuntimeError:
                print "Warning: adjust termination_probability"
                smi = "CC(CC(CCC))(CC)CC"
            self.assertOne(smi)

    def test_bonds(self):
        "Test that all four types of basic bond work in a simple string."
        smi = "CCC%sCCCC"
        smileses = [smi % bondsym for bondsym in BASIC_BONDS.keys()]
        self.assertMany(smileses)

    def test_bonds_tiny(self):
        "Test that all four types of basic bond work in a tiny string."
        smi = "C%sC"
        smileses = [smi % bondsym for bondsym in BASIC_BONDS.keys()]
        self.assertMany(smileses)

    def test_dots(self):
        "Test that dots work correctly."
        smileses = [
            r"CC(.CC)CCCC",
            r"CC(C.CC)CCCC",
            r"C(C(C(C(C.CCCC)C)C)C)CC",
            r"C(C(C(C(.CCCC)C)C)C)CC",
            r"CCCC.CCCCC",
        ]
        self.assertMany(smileses)

    def test_rings(self):
        "Test that basic usage of ringlinks is functional."
        smileses = [
            r"C1CCCC1",
            r"CCC2CC(CC(CC)(C2)CC)CC",
            r"CCC(CC)(CC1)CCCC1",
            r"C1(C1)",
        ]
        self.assertMany(smileses)

    def test_rings_many(self):
        "Test that multiple ringlinks can be used at the same time."
        smileses = [
            r"CCCC1CC2CCCC1C2",
            r"CCC1CC2CCCC2C1",
            r"CCC1CC2CC3CC1CC2CC3",
            r"CC2CCCCC3CCCC1CC2C3C1CC",
            r"C1C2C3C4C5C6C7C8C9C%10C%11CCC%11C%10C9C8C7C6C5C4C3C2C1",
            r"C123CCC1C2C3",
            r"CC2(CC(CC2)C2CCC2C)C",
            r"C2(CCC2CCC2CCC2)",
        ]
        self.assertMany(smileses)

    def test_rings_percent(self):
        "Test that percent-escaping ringlinks works."
        smileses = [
            r"CCC%10CCCC%10",
            r"CCC%99CCC%22CCC%22CCCC%99",
            r"CCC1CC%12CCC1CC%12CCC",
            r"CCC1CC=%12CCC=1CC%12CCC",
        ]
        self.assertMany(smileses)

    def test_rings_with_bonds(self):
        "Test that prefacing ringlinks with bonds works."
        smileses = [
            r"CCC=1CCC1",
            r"CCC=1CCCC=1",
            r"C1CCCC#2CCC1C2C",
            r"C1CCCC#%88CCC1C%88C",
            r"C1CCCC=%88CCC1C=%88C",
        ]
        self.assertMany(smileses)

    def test_dotteds(self):
        "Test that dots work correctly."
        smiles_groups = [
            (r"CCC.CCCC", [r"CCC", r"CCCC"]),
            (r"C.CC", [r"C", r"CC"]),
            (r"C.CC.CCC", [r"C", r"CC", r"CCC"]),
            (r"C(CC.CCC)CCC.CCCC", [r"CCC", r"C(CC)CCC", r"CCCC"])
        ]
        for (dotted, split) in smiles_groups:
            first = to_canonical('.'.join([to_canonical(i) for i in split]))
            second = smilesify(moleculify(dotted))
            self.assertEqual(first, second)

    def test_all_carbon(self):
        "Summatively test that rings, bonds, and branching all work."
        smileses = [
            r"C1C(C(C(C)(C)C)C=C(C#C)C2)C2C(C)(C)C=1",
        ]
        self.assertMany(smileses)

    def test_examples_easy(self):
        "Test all the advanced example smiles."
        smileses = example_smiles_easy
        self.assertMany(smileses)

    @unittest.skip("Not implemented yet")
    def test_examples_cistrans(self):
        "Test all the advanced example smiles."
        smileses = example_smiles_cistrans
        self.assertMany(smileses)

    def test_examples_chiral(self):
        "Test all the advanced example smiles."
        smileses = example_smiles_chiral
        self.assertMany(smileses)

    @unittest.skip("Not implemented yet")
    def test_examples_hard(self):
        "Test all the advanced example smiles."
        smileses = example_smiles_hard
        self.assertMany(smileses)

    def test_examples_hard_no_chirality_or_cistrans(self):
        "Test all the advanced example smiles."
        smileses = example_smiles_hard_no_chirality_or_cistrans
        self.assertMany(smileses)

    def test_examples_hard_no_cistrans(self):
        "Test all the advanced example smiles."
        smileses = example_smiles_hard_no_cistrans
        self.assertMany(smileses)

    def test_examples_atoms(self):
        "Test all the advanced example smiles."
        smileses = example_smiles_atoms
        self.assertMany(smileses)

    def test_wacky_hydrogens(self):
        "Test that implicit/explicit hydrogens work correctly."
        self.assertSameMany([
            r'C',
            r'[CH4]',
            r'[H][C]([H])([H])[H]',
            r'[C]([H])([H])([H])[H]',
        ])
        self.assertSameMany([
            r'[C][H]',
            r'[CH]',
            r'[CH1]',
            r'[H][C]',
        ])
        self.assertSameMany([
            r'[C]([H])[H]',
            r'[H][C][H]',
            r'[CH2]',
        ])
        self.assertSameMany([
            r'[C]([H])([H])[H]',
            r'[H][C]([H])[H]',
            r'[CH3]',
        ])
        self.assertSameMany([
            r'[C]([H])([H])([H])([H])[H]',
            r'[H][C]([H])([H])([H])[H]',
            r'[CH5]',
        ])
        self.assertSameMany([
            r'[H][C]([H])([H])([H])([H])[H]',
            r'[CH6]',
        ])
        self.assertMany([
            r'[C]',
            r'[CH7]',
            r'[CH8]',
            r'[CH9]',
        ])
        self.assertSameMany([
            r'[SeH]O',
            r'O[Se][H]',
        ])
        self.assertSameMany([
            r'[ClH]N',
            r'N[Cl][H]',
        ])
        self.assertSameMany([
            r'[Br][F]',
            r'FBr',
        ])
        self.assertSameMany([
            r'PI',
            r'P[I]',
        ])
        self.assertSameMany([
            r'OO',
            r'[OH]O',
            r'[OH1][OH]',
            r'[H][O][O][H]',
            r'[O]([H])[O][H]',
        ])

    def test_chirality_equivalence(self):
        """
        Test that chirality works PERFECTLY.
        These are only a few of the many ways to represent this chiral molecule.
        """
        self.assertSameMany([
            r'Br[C@](O)(N)C',
            r'Br[C@](C)(O)N',
            r'Br[C@](N)(C)O',
            r'N[C@](Br)(O)C',
            r'O[C@](Br)(C)N',
            r'C[C@](Br)(N)O',
            r'[C@](N)(Br)(O)C',
            r'[C@](N)(O)(C)Br',
            r'[C@](N)(C)(Br)O',
            r'[C@](O)(Br)(C)N',
            r'[C@](O)(C)(N)Br',
            r'[C@](O)(N)(Br)C',
            r'[C@](Br)(O)(N)C',
            r'[C@](Br)(C)(O)N',
            r'[C@](Br)(N)(C)O',
            r'[C@](C)(Br)(N)O',
            r'[C@](C)(N)(O)Br',
            r'[C@](C)(O)(Br)N',
            r'C[C@@](Br)(O)N',
            r'Br[C@@](N)(O)C',
            r'[C@@](C)(Br)(O)N',
            r'[C@@](Br)(N)(O)C',
        ])


    def test_ethanol(self):
        "All of these are equivalent to ethanol, C2H5OH."
        self.assertSameMany([
            r'CCO'
            r'OCC'
            r'C(O)C'
            r'[CH3][CH2][OH]'
            r'[H][C]([H])([H])C([H])([H])[O][H]'
        ])

    def test_benzene(self):
        "All of these are equivalent to benzene, C6H6."
        self.assertSameMany([
            r'c1ccccc1',
            r'c:1:c:c:c:c:c1',
            r'c1:c:c:c:c:c:1',
            r'c1:cc:c:cc:1',
            r'C:1:C:C:C:C:C1',
            r'C1:C:C:C:C:C:1',
            r'C:1:C:C:C:C:C:1',
            r'C1=CC=CC=C1',
            r'C=1C=CC=CC1',
            r'C=1C=CC=CC=1',
        ])

    def test_phenol(self):
        "All of these are equivalent to phenol, C6H5OH."
        self.assertSameMany([
            r'Oc1ccccc1',
            r'c1ccccc1O',
            r'c1(O)ccccc1',
            r'c1(ccccc1)O',
        ])

    def test_aromaticity(self):
        "Some diverse tests of aromaticity notation."
        self.assertMany([
            r'c1ccccc1-c2ccccc2',
            r'c1cc(F)c(CC)cc1',
        ])
        self.assertSame(r'C1=CC=CC(CCC2)=C12', r'c1ccc2CCCc2c1')
        self.assertSame(r'c1occc1', r'C1OC=CC=1')

        ## cyclobutadiene does not canonicalize correctly
        ## https://sourceforge.net/p/openbabel/bugs/940/
        # self.assertSameMany([
        #     r'c1ccc1',
        #     r'C1=CC=C1',
        # ])

    def test_truly_evil_molecule(self):
        "Why must I hold my code to these standards..."
        self.assertOne("O1CCN[C@]11(CCCO1)")

    def test_chirality_frustration(self):
        self.assertOne("[C@](Br)(F)(O)[H]")

    # def test_nonstandard(self):
    #     "Nonstandard forms of SMILES that we could encounter."
    #     self.assertSame("C((C))O", "C(C)O")
    #     self.assertSame("(N1CCCC1)", "N1CCCCC1")
    #     self.assertSame("[Na+]..[Cl-]", "[Na+].[Cl-]")
    #     self.assertSame(".CCO", "CCO")
    #     self.assertSame("CCO.", "CCO")
    #     self.assertSame("C1CCC", "CCCC")
    #     self.assertSame("D[CH3]", "[2H][CH3]")
    #     self.assertSame("T[CH3]", "[3H][CH3]")
