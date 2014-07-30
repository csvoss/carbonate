"""
toMolecule.py

This code contains the method moleculify(smiles), which converts a SMILES string to a Molecule object:
http://en.wikipedia.org/wiki/SMILES

This code is not yet complete.

Only ``moleculify'' is meant to be public-facing. The other methods in this file are private.

This would allow us to take in SMILES strings as inputs, from students or professors, to specify their own molecules to use in problems. Hooray!

Refer to:
    https://github.com/alex/rply
    https://pypi.python.org/pypi/rply/0.5.1
    http://www.opensmiles.org/opensmiles.html
"""

from molecularStructure import *
import string
from rply import ParserGenerator, LexerGenerator, ParsingError
from rply.token import BaseBox
from toSmiles import smilesify, to_canonical
import unittest
from random import random, randrange


debug = True

##### HELPER FUNCTIONS #####

BASIC_BONDS = {'-':1, '=':2, '#':3, '$':4}

# def plus_it(name, classname):
#     """
#     Create some productions representing <name>+.
#     The productions output a list of <classname> objects.
#     name :: str. example: "ringbond"
#     classname :: type. example: Ringbond
#     plusname = name + "plus"
#     """
#     plusname = name + "plus"
#     @pg.production("%s : %s" % (plusname, name))
#     def plus_empty(p):
#         return []
#     @pg.production("%s : %s %s" % (plusname, plusname, name))
#     def plus_full(p):
#         xplus = p[0]
#         x = p[1]
#         assert_isinstance(xplus, list)
#         assert_isinstance(x, classname)
#         return xplus + [x]
#     return plus_empty, plus_full

def bond_at(bond, m1, a1, a2, m2):
    """
    Create a bond between molecules m1 and m2, at atoms a1 and a2, using bond.
    THIS METHOD IS INCOMPLETE. #TODO
    """
    if bond == '.':
        raise StandardError("bond_at called with '.' as bond")
    elif bond == 'default':
        ## TODO: support aromaticity!!!
        m1.addMolecule(m2, a2, a1, 1)
    elif bond in BASIC_BONDS.keys():
        ## connect them together
        # raise StandardError('=') ## TODO
        m1.addMolecule(m2, a2, a1, BASIC_BONDS[bond])
    else:
        raise NotImplementedError(":, /, and \\ are not-yet-implemented bond types.")
        ## TODO.
        ## Idea: Implement / and \\ by adding them as relevant flags,
        ## then postprocessing molecules to convert them into cis/trans centers

def debug_decorator(f):
    """
    Apply this decorator to a function to make that function print at you
    whenever it is called.
    """
    def new_f(*args, **kwargs):
        print "----- Calling %s with %s and %s" % (f, args, kwargs)
        f(*args, **kwargs)
    return new_f


def assert_isinstance(i, t):
    """
    Raise an error if type(i) is not t.
    """
    assert isinstance(i, t), "In parser: Object %s is of type %s; should be %s" % (repr(i), repr(type(i)), repr(t))

##### LEXER #####

## Create a Lexer, which splits up a string into tokens.
## Tokens may be, for example, SYMBOLs or LETTERs or TERMINATORs or DIGITs.
## Whitespace will be ignored.

lg = LexerGenerator()
lg.ignore(r"\s+")

#lg.add('SYMBOL', r'[@#$%\*\(\)\[\]=\+\-:/\\\.]') ## SPLIT OUT this one
SYMBOLS_DICT = {'@':'@','%':'%','*':'\*','(':'\(',')':'\)','[':'\[',']':'\]','+':'\+','.':'\.'}
SYMBOLS = SYMBOLS_DICT.keys()
for sym in SYMBOLS:
    lg.add(sym, SYMBOLS_DICT[sym])

#lg.add('LETTER', r'[A-IK-PRSTXYZa-ik-pr-vy]') ## SPLIT OUT this one
LETTERS = ["C"] ## TODO this is for debugging
#LETTERS = [c for c in "ABCDEFGHIKLMNOPRSTXYZabcdefghiklmnoprstuvy"] # omissions intentional
for sym in LETTERS:
    lg.add(sym, sym)

lg.add('DIGIT', r'[0-9]')
lg.add('TERMINATOR', r'[ \t\r\n]')
lg.add('BOND', r'[\-=#\$:/\\]')




##### PARSER #####
# www.opensmiles.org/opensmiles.html

pg = ParserGenerator(
    SYMBOLS + LETTERS + ['DIGIT', 'TERMINATOR', 'BOND'],
    precedence = [
    ],
    cache_id='molparser',
)

# main :: [Molecule].
@pg.production("main : smiles")
def main_production(p):
    return p[0]

# smiles ::= terminator | chain terminator
# smiles :: [Molecule].
@pg.production("smiles : chain terminator")
@pg.production("smiles : chain")
def smiles_production(p):
    chain = p[0]
    assert_isinstance(chain, Chain)
    ## TODO: Post-processing of ringbond list at chain.ring_data_list
    ## TODO: Post-processing of aromaticity
    # print smilesify(chain.dotted_molecules)
    return [chain.molecule] + chain.dotted_molecules
@pg.production("smiles : terminator")
def smiles_empty(p):
    return []

# terminator ::= SPACE TAB | LINEFEED | CARRIAGE_RETURN | END_OF_STRING
# terminator :: None.
@pg.production("terminator : TERMINATOR")
def terminator_production(p):
    return None


##### BONDS #####
# http://www.opensmiles.org/opensmiles.html#bonds

class Chain(BaseBox):
    def __init__(self, m, first, last):
        """m :: Molecule.
        first, last :: Atom."""
        self.molecule = m
        self.first_atom = first
        self.last_atom = last
        self.dotted_molecules = []
        self.ring_data_list = []

    def append_dotteds(self, molecules):
        "molecule :: [Molecule]."
        for molecule in molecules:
            self.dotted_molecules.append(molecule)

    def append_ring(self, ring_index, bond_char):
        """ring_index :: int.
        bond_char :: str."""
        self.ring_data_list.append((self, ring_index, bond_char))

    def append_rings(self, ring_list):
        """ring_list :: [(Atom, int ring_index, str bond_char)]."""
        self.ring_data_list += ring_list


# chain ::= branched_atom | chain branched_atom | chain bond branched_atom | chain dot branched_atom
# chain :: Chain.
@pg.production("chain : branched_atom")
def chain_production(p):
    m, a = p[0].as_tuple()
    assert_isinstance(m, Molecule)
    assert_isinstance(a, Atom)
    return Chain(m, a, a)
@pg.production("chain : chain unparenbranch")
def chain_unparenbranch_testing_temporary(p): #TODO: Does it work?
    chain = p[0]
    branch = p[1]
    assert_isinstance(chain, Chain)
    assert_isinstance(branch, Branch)
    ## Now, ligate branch onto chain at the last atom of chain
    bond_to_use = branch.bond_or_dot
    is_dot = (bond_to_use == '.')
    m1, a1 = chain.molecule, chain.last_atom
    m2, a2 = branch.chain.molecule, branch.chain.first_atom
    assert_isinstance(m1, Molecule)
    assert_isinstance(a1, Atom)
    assert_isinstance(m2, Molecule)
    assert_isinstance(a2, Atom)
    if is_dot:
        chain.append_dotteds([m2])
        chain.append_dotteds(branch.chain.dotted_molecules)
        chain.append_rings(branch.chain.ring_data_list) ## Atom, int, str -> int, char
        return chain
    else:
        bond_at(bond_to_use, m1, a1, a2, m2)
        chain.last_atom = branch.chain.last_atom
        chain.append_dotteds(branch.chain.dotted_molecules)
        chain.append_rings(branch.chain.ring_data_list)
        return chain

# @pg.production("chain : chain branched_atom")
# def chain_chain(p):
#     return chain_bond([p[0], 'default', p[1]])
# @pg.production("chain : chain bond branched_atom")
# def chain_bond(p):
#     bond = p[1]
#     assert_isinstance(bond, basestring)
#     chain = p[0]
#     m1, a1 = chain.molecule, chain.last_atom
#     m2, a2 = p[2].as_tuple()
#     assert_isinstance(m1, Molecule)
#     assert_isinstance(a1, Atom)
#     assert_isinstance(m2, Molecule)
#     assert_isinstance(a2, Atom)
#     bond_at(bond, m1, a1, a2, m2)
#     chain.last_atom = a2
#     chain.append_dotteds(p[2].dotted_molecules)
#     chain.append_rings(p[2].ring_data_list)
#     return chain
# @pg.production("chain : chain dot branched_atom")
# def chain_dot(p):
#     chain = p[0]
#     m1, a1 = chain.molecule, chain.last_atom
#     dot = p[1]
#     m2, a2 = p[2].as_tuple()
#     assert dot == '.'
#     assert_isinstance(m1, Molecule)
#     assert_isinstance(a1, Atom)
#     assert_isinstance(m2, Molecule)
#     assert_isinstance(a2, Atom)
#     chain.append_dotteds([m2])
#     chain.append_dotteds(p[2].dotted_molecules)
#     chain.append_rings(p[2].ring_data_list) ## Atom, int, str -> int, char
#     return chain


# bond ::= '-' | '=' | '#' | '$' | ':' | '/' | '\\'
# bond :: str
# @pg.production("bond : -")
# @pg.production("bond : =")
# @pg.production("bond : #")
# @pg.production("bond : $")
# @pg.production("bond : colon")
# @pg.production("bond : /")
# @pg.production("bond : \\")
# def bond_production(p):
#     bond = p[0].getstr()
#     assert_isinstance(bond, basestring)
#     return bond
@pg.production("bond : BOND")
def bond_production(p):
    bond = p[0].getstr()
    assert_isinstance(bond, basestring)
    return bond

# dot ::= '.'
# dot :: str
@pg.production("dot : .")
def dot_production(p):
    dot = p[0].getstr()
    assert_isinstance(dot, basestring)
    return dot

class BranchedAtom(BaseBox):
    def __init__(self, m, a):
        """m :: Molecule.
        a :: Atom."""
        self.molecule = m
        self.atom = a
        self.dotted_molecules = []
        self.ring_data_list = [] ## :: [(Atom, int ring_index, str bond_char)].

    def as_tuple(self):
        "return :: (Molecule, Atom)."
        return (self.molecule, self.atom)

    def append_dotteds(self, molecules):
        "molecule :: [Molecule]"
        for molecule in molecules:
            self.dotted_molecules.append(molecule)

    def append_ring(self, ring_index, bond_char):
        """ring_index :: int.
        bond_char :: str."""
        self.ring_data_list.append((self, ring_index, bond_char))

    def append_rings(self, ring_list):
        """ring_list :: [(Atom, int ring_index, str bond_char)]."""
        self.ring_data_list += ring_list


def add_rings_and_branches(atom, ringbonds, branches):
    "return :: BranchedAtom"
    assert_isinstance(atom, Atom)
    assert_isinstance(ringbonds, list)
    assert_isinstance(branches, list)
    output = BranchedAtom(Molecule(atom), atom)
    for branch in branches:
        output.append_dotteds(branch.chain.dotted_molecules)
        if branch.bond_or_dot == '.':
            pass
        else:
            other_molecule = branch.chain.molecule
            other_atom = branch.chain.first_atom
            bond_at(
                branch.bond_or_dot,
                output.molecule,
                output.atom,
                other_atom,
                other_molecule
            )
    for ringbond in ringbonds:
        output.append_ring(ringbond.index, ringbond.bond)
    return output

# branched_atom :: BranchedAtom.
# branched_atom2 :: Atom, [Ringbond], [Branch].
# ringed_atom :: Atom, [Ringbond], [Branch].
@pg.production("branched_atom : branched_atom2")
def branched_to_branched2(p):
    " :: BranchedAtom"
    atom, ringbonds, branches = p[0]
    return add_rings_and_branches(atom, ringbonds, branches)
# @pg.production("branched_atom : branched_atom2 unparenbranch")  ## TODO I think this rule is unnecessary.
# def branched_to_branched2_with_unparenbranch(p):
#     " :: BranchedAtom"
#     atom, ringbonds, branches = p[0]
#     unparenbranch = p[1]
#     branches = branches + [unparenbranch]
#     return add_rings_and_branches(atom, ringbonds, branches)
@pg.production("branched_atom2 : ringed_atom")
def branched_to_ringed(p):
    " :: Atom, [Ringbond], [Branch]"
    return p[0]
# @pg.production("branched_atom2 : atom")
# def branched_to_plain_atom_maybe_this_will_work(p):
#     " :: Atom, [Ringbond], [Branch]"
#     atom = p[0]
#     return atom, [], []
@pg.production("branched_atom2 : branched_atom2 branch")
def branched_to_branched(p):
    " :: Atom, [Ringbond], [Branch]"
    atom, ringbonds, branches = p[0]
    return atom, ringbonds, branches+[p[1]]
@pg.production("ringed_atom : ringed_atom ringbond")
def ringed_to_ringed(p):
    " :: Atom, [Ringbond], [Branch]"
    atom, ringbonds, branches = p[0]
    return atom, ringbonds+[p[1]], branches
@pg.production("ringed_atom : atom")
def ringed_to_atom(p):
    " :: Atom, [Ringbond], [Branch]"
    atom = p[0]
    return atom, [], []


class Ringbond(BaseBox):
    def __init__(self, index, bond):
        assert_isinstance(index, int)
        assert_isinstance(bond, basestring)
        self.index = index
        self.bond = bond
class Branch(BaseBox):
    def __init__(self, bond_or_dot, chain):
        self.bond_or_dot = bond_or_dot
        self.chain = chain

# branch ::= '(' chain ')' | '(' bond chain ')' | '(' dot chain ')'
# branch :: Branch.
# unparenbranch :: Branch.
@pg.production("branch : ( unparenbranch )")
def unbranch_production(p):
    assert p[0].getstr()=='(' and p[2].getstr()==')', \
        "Not enclosed in parentheses: %s" % str(p)
    # print smilesify(p[1].chain.molecule)
    return p[1]
@pg.production("unparenbranch : chain")
def branch_production(p):
    chain = p[0]
    return Branch('-', chain)
@pg.production("unparenbranch : bond chain")
def branch_bond(p):
    bond = p[0]
    chain = p[1]
    assert_isinstance(bond, basestring)
    assert_isinstance(chain, Chain)
    return Branch(bond, chain)
@pg.production("unparenbranch : dot chain")
def branch_dot(p):
    dot = p[0]
    chain = p[1]
    assert_isinstance(dot, basestring)
    assert_isinstance(chain, Chain)
    return Branch(dot, chain)

# ringbond ::= bond? DIGIT | bond? '%' DIGIT DIGIT
# ringbond :: Ringbond
@pg.production("ringbond : ringbond2")  ## TODO aaaughh
def ringbond_to_ringbond2(p):
    return p[0]
@pg.production("ringbond : DIGIT")
def ringbond_prod_1(p):
    return ringbond_prod_2(['default'] + p)
@pg.production("ringbond2 : bond DIGIT")
def ringbond_prod_2(p):
    bond = p[0]
    index = int(p[1].getstr())
    return Ringbond(index, bond)
@pg.production("ringbond : percentdigit")
def ringbond_prod_3(p):
    return Ringbond(p[0], 'default')
# @pg.production("ringbond2 : bond percentdigit")  ## TODO uncomment
# def ringbond_prod_4(p):
#     bond = p[0]
#     return Ringbond(p[1], bond)

# percentdigit : int
@pg.production("percentdigit : % DIGIT DIGIT")
def percentdigit(p):
    useless_percent_sign = p[0]
    index_1 = int(p[1].getstr())
    index_2 = int(p[2].getstr())
    index = 10*index_1 + index_2
    return index

##### ATOMS #####

# atom ::= bracket_atom | aliphatic_organic | aromatic_organic | '*'
@pg.production("atom : C")
def atom_production(p):
    ## TODO -- Temporary
    return Atom("C")

##### ORGANIC SUBSET ATOMS #####

# aliphatic_organic ::= 'B' | 'C' | 'N' | 'O' | 'S' | 'P' | 'F' | 'Cl' | 'Br' | 'I'

# aromatic_organic ::= 'b' | 'c' | 'n' | 'o' | 's' | 'p'

##### BRACKET ATOMS #####

# bracket_atom ::= '[' isotope? symbol chiral? hcount? charge? class? ']'

# symbol ::= element_symbols | aromatic_symbols | '*'

# isotope ::= NUMBER

# element_symbols ::= 'H' | 'He' | 'Li' | 'Be' | 'B' | 'C' | 'N' | 'O' | 'F' | 'Ne' | 'Na' | 'Mg' | 'Al' | 'Si' | 'P' | 'S' | 'Cl' | 'Ar' | 'K' | 'Ca' | 'Sc' | 'Ti' | 'V' | 'Cr' | 'Mn' | 'Fe' | 'Co' | 'Ni' | 'Cu' | 'Zn' | 'Ga' | 'Ge' | 'As' | 'Se' | 'Br' | 'Kr' | 'Rb' | 'Sr' | 'Y' | 'Zr' | 'Nb' | 'Mo' | 'Tc' | 'Ru' | 'Rh' | 'Pd' | 'Ag' | 'Cd' | 'In' | 'Sn' | 'Sb' | 'Te' | 'I' | 'Xe' | 'Cs' | 'Ba' | 'Hf' | 'Ta' | 'W' | 'Re' | 'Os' | 'Ir' | 'Pt' | 'Au' | 'Hg' | 'Tl' | 'Pb' | 'Bi' | 'Po' | 'At' | 'Rn' | 'Fr' | 'Ra' | 'Rf' | 'Db' | 'Sg' | 'Bh' | 'Hs' | 'Mt' | 'Ds' | 'Rg' | 'Cn' | 'Fl' | 'Lv' | 'La' | 'Ce' | 'Pr' | 'Nd' | 'Pm' | 'Sm' | 'Eu' | 'Gd' | 'Tb' | 'Dy' | 'Ho' | 'Er' | 'Tm' | 'Yb' | 'Lu' | 'Ac' | 'Th' | 'Pa' | 'U' | 'Np' | 'Pu' | 'Am' | 'Cm' | 'Bk' | 'Cf' | 'Es' | 'Fm' | 'Md' | 'No' | 'Lr'

# aromatic_symbols ::= 'b' | 'c' | 'n' | 'o' | 'p' | 's' | 'se' | 'as'

##### CHIRALITY #####

# chiral ::= '@' | '@@' | '@TH1' | '@TH2' | '@AL1' | '@AL2' | '@SP1' | '@SP2' | '@SP3' | '@TB1' | '@TB2' | '@TB3' | ... | '@TB20' | '@OH1' | '@OH2' | '@OH3' | ... | '@OH30' | '@TB' DIGIT DIGIT | '@OH' DIGIT DIGIT

##### HYDROGENS #####

# hcount ::= 'H' | 'H' DIGIT

##### CHARGES #####

# charge ::= '-' | '-' DIGIT? DIGIT | '+' | '+' DIGIT? DIGIT | '--' deprecated | '++' deprecated

##### ATOM CLASS #####

# class ::= ':' NUMBER



def dictof(thing):
    try:
        return str(thing.__dict__)
    except:
        return str(thing)

@pg.error
def error_handler(token, expected=[]):
    raise ValueError("Ran into a %s (%s) where it wasn't expected. %s. Instead expected: %s." % (repr(token.name), repr(token.value), dictof(token.source_pos), repr(expected)))



lexer = lg.build()
parser = pg.build()






def moleculify(smiles):
    """
    smiles :: str or [str]. SMILES string(s) e.g. "CC(CN)CCC(O)O"
    return :: [Molecule].

    Raises a StandardError if the SMILES string contains
    as-yet-unsupported features (like delocalization).
    """
    #return example_molecule()
    if isinstance(smiles, list):
        return map(moleculify, smiles)
    else:
        try:
            lexed = lexer.lex(smiles)
            return parser.parse(lexed)
        except ParsingError as e:
            raise StandardError(e.getsourcepos())


def example_molecule():
    """
    This method intended as a code example of how to create and return a molecule.
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
    mol4.addBond(c42,c43,1)

    return mol4.withHydrogens()

def example_smiles():
    """
    return :: [SMILES str].
    """
    return [
        "C(/F)(\Cl)=C(/Br)\I",  ## cis/trans test
        "F\C=C(/Br)\I",     ## more cis/trans
        "F\C=C/Br",         ## more cis/trans
        "C1CCCCC1",         ## cyclohexane
        "F/C=C/F",          ## trans-difluoroethene
        "F/C=C\F",          ## cis-difluoroethene
        "N[C@@H](C)C(=O)O", ## L-alanine
        "N[C@H](C(=O)O)C",  ## also L-alanine
        "N[C@H](C)C(=O)O",  ## D-alanine (less common)
        "[2H]C(Cl)(Cl)Cl",  ## deuterochloroform (hydrogen-2)
        "N#N",              ## dinitrogen
        "CN=C=O",           ## methyl isocyanate
        "[Cu+2].[O-]S(=O)(=O)[O-]" ## copper(ii) sulfate
        "CCC[C@@H](O)CC\C=C\C=C\C#CC#C\C=C\CO", ## enanthotoxin
        "COC(=O)C(\C)=C\C1C(C)(C)[C@H]1C(=O)O[C@@H]2C(C)=C(C(=O)C2)CC=CC=C", ## pyrethrin II
        "OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H](O)[C@@H](O)1", ## glucose
        "CC(=O)OCCC(/C)=C\C[C@H](C(C)=C)CCC=C", ## some pheromone or another
        "CC[C@H](O1)CC[C@@]12CCCO2", ## another pheromone
        "CC(C)[C@@]12C[C@@H]1[C@@H](C)C(=O)C2", ## alpha-thujone
        "N1CCN(CC1)C(C(F)=C2)=CC(=C2C4=O)N(C3CC3)C=C4C(=O)O",
        "C[C@@](C)(O1)C[C@@H](O)[C@@]1(O2)[C@@H](C)[C@@H]3CC=C4[C@]3(C2)C(=O)C[C@H]5[C@H]4CC[C@@H](C6)[C@]5(C)CC(N7)C6NC(C[C@@]89(C))C7C[C@@H]8CC[C@@H]%10[C@@H]9C[C@@H](O)[C@@]%11(C)C%10=C[C@H](O%12)[C@]%11(O)[C@H](C)[C@]%12(O%13)[C@H](O)C[C@@]%13(C)CO", ## Cephalostatin-1 without aromaticity. Note the % in front of ring closure labels above 9.
    ]



##### UNIT TESTS #####

def two_smiles(smiles):
    ret = [to_canonical(smiles)], smilesify(moleculify(smiles))
    # print ret
    return ret

class TestBasic(unittest.TestCase):
    def test_carbons(self):
        for num in range(1, 21):
            smi = "C" * num
            self.assertEqual(*two_smiles(smi))

class TestBranches(unittest.TestCase):
    def test_branches(self):
        termination_probability = 0.90
        max_run_length = 5
        def random_branching_carbons():
            def end_or_branch():
                if random() < termination_probability:
                    return "C" * (1 + randrange(max_run_length))
                else:
                    return random_branching_carbons()
            a = end_or_branch()
            b = end_or_branch()
            c = end_or_branch()
            d = end_or_branch()
            return a + "(" + b + ")(" + c + ")" + d

        for num in range(1,20):
            try:
                smi = random_branching_carbons()
            except RuntimeError:
                print "Warning: adjust termination_probability"
                smi = "CC(CC(CCC))(CC)CC"
            # print smi ##TODO remove
            self.assertEqual(*two_smiles(smi))

class TestBonds(unittest.TestCase):
    def test_bonds(self):
        smi = "CCC%sCCCC"
        for bondsym in BASIC_BONDS:
            self.assertEqual(*two_smiles(smi % bondsym))

    def test_bonds2(self):
        smi = "C%sC"
        for bondsym in BASIC_BONDS:
            self.assertEqual(*two_smiles(smi % bondsym))

class TestRings(unittest.TestCase):
    def test_rings(self):
        smi = "CCCCCC"
        self.assertEqual(*two_smiles(smi))



## TODO temporary
m = moleculify
s = smilesify