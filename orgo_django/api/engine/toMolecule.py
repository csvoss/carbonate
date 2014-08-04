"""
toMolecule.py

This code contains the method moleculify(smiles), which converts a SMILES
string to a Molecule object:
http://en.wikipedia.org/wiki/SMILES

This code is not yet complete.

Only ``moleculify'' is meant to be public-facing. The other methods in
this file are private.

This would allow us to take in SMILES strings as inputs, from students or
professors, to specify their own molecules to use in problems. Hooray!

Refer to:
    https://github.com/alex/rply
    https://pypi.python.org/pypi/rply/0.5.1
    http://www.opensmiles.org/opensmiles.html
"""
from random import random, randrange
import itertools
import re
import unittest

from rply import ParserGenerator, LexerGenerator, ParsingError
from rply.token import BaseBox

from molecularStructure import Molecule, Atom
from toSmiles import smilesify, to_canonical

DEBUG = True

############################
##### PUBLIC FUNCTIONS #####
############################

def moleculify(smiles):
    """
    smiles :: str or [str]. SMILES string(s) e.g. "CC(CN)CCC(O)O"
    return :: [Molecule].

    Raises a StandardError if the SMILES string contains
    as-yet-unsupported features (like delocalization).
    """
    #return example_molecule()
    if isinstance(smiles, list):
        return [moleculify(i) for i in smiles]
    else:
        try:
            smiles = preprocess(smiles)
            lexed = LEXER.lex(smiles)
            return PARSER.parse(lexed)
        except ParsingError as e:
            raise StandardError(e.getsourcepos())


############################
##### HELPER FUNCTIONS #####
############################

BASIC_BONDS = {'-':1, '=':2, '#':3, '$':4}

DEFAULT_BOND = 'default'

def assert_isinstance(instance, some_type):
    """
    Raise an error if type(instance) is not some_type.
    """
    assert isinstance(instance, some_type), \
        "In parser: Object %s is of type %s; should be %s" % \
        (repr(instance), repr(type(instance)), repr(some_type))

def bond_at(bond, mol1, atom1, atom2, mol2):
    """
    Create a bond between molecules mol1 and mol2, at atoms atom1 and atom2,
    using bond.

    REQUIRES that atom1 precede atom2 in the SMILES string.

    mol1, mol2 :: Molecule.
    atom1, atom2 :: Atom.
    return :: None. Works via side effects.
    """
    if bond == '.':
        raise StandardError("bond_at called with '.' as bond")
    elif bond == DEFAULT_BOND:
        ## TODO: support aromaticity!!!
        mol1.addMolecule(mol2, atom2, atom1, 1)
    elif bond in BASIC_BONDS.keys():
        ## connect them together
        mol1.addMolecule(mol2, atom2, atom1, BASIC_BONDS[bond])
    else:
        raise NotImplementedError(":, /, and \\ are not yet implemented.")
        ## TODO.
        ## Idea: Implement / and \\ by adding them as relevant flags,
        ## then postprocessing molecules to add cis/trans centers

def dictof(thing):
    "Printable `thing.__dict__` if possible, else just `thing`"
    try:
        return str(thing.__dict__)
    except AttributeError:
        return str(thing)

def debug_decorator(func):
    """
    Apply this decorator to a function to make that function print at you
    whenever it is called.
    """
    def new_f(*args, **kwargs):
        "The "
        print "----- Calling %s with %s and %s" % (func, args, kwargs)
        func(*args, **kwargs)
    return new_f

def int_from_digit(digit):
    """
    digit :: Token.
    return :: int.
    """
    return int(digit.getstr())


#################
##### LEXER #####
#################

## Create a Lexer, which splits up a string into tokens.
## Tokens may be, for example, SYMBOLs or LETTERs or TERMINATORs or DIGITs
## Whitespace will be ignored.

LG = LexerGenerator()
LG.ignore(r"\s+")

SYMBOLS_DICT = {
    '@': '@', '%': '%', '*': r'\*', '(': r'\(', ')': r'\)',
    '[': r'\[', ']': r'\]', '+': r'\+', '.': r'\.'
}
SYMBOLS = SYMBOLS_DICT.keys()
for sym in SYMBOLS:
    LG.add(sym, SYMBOLS_DICT[sym])

LETTERS = [letter for letter in "ABCDEFGHIKLMNOPRSTUVWXYZabcdefghiklmnoprstuvy"]
for sym in LETTERS:
    LG.add(sym, sym)

BOND_SYMBOLS_DICT = {
    '-': r'\-',
    '=': r'=',
    '#': r'#',
    '$': r'\$',
    'colon': r':',
    '/': r'/',
    '\\': r'\\',
}
for sym in BOND_SYMBOLS_DICT:
    LG.add(sym, BOND_SYMBOLS_DICT[sym])
BOND_SYMBOLS = BOND_SYMBOLS_DICT.keys()

BOND_SYMBOLS_TILDE_DICT = {
    '~-~': r'~\-~',
    '~=~': r'~=~',
    '~#~': r'~#~',
    '~$~': r'~\$~',
    '~colon~': r'~:~',
    '~/~': r'~/~',
    '~\\~': r'~\\~',
}
for sym in BOND_SYMBOLS_TILDE_DICT:
    LG.add(sym, BOND_SYMBOLS_TILDE_DICT[sym])
BOND_SYMBOLS_TILDE = BOND_SYMBOLS_TILDE_DICT.keys()


LG.add('DIGIT', r'[0-9]')
LG.add('TERMINATOR', r'[ \t\r\n]')
# LG.add('BOND', r'[\-=#\$:/\\]')

FOR_PREPROCESSOR = {
    r'\-': '~-~',
    r'=': '~=~',
    r'#': '~#~',
    r'\$': '~$~',
    r':': '~:~',
    r'/': '~/~',
    r'\\': '~\\~',
}
ATOM_MATCHER = r'\[ABCDEFGHIKLMNOPRSTUVWXYZabcdefghiklmnoprstuvy'

def preprocess(smiles):
    """
    Okay, here's the deal. We can't have double bonds to rings and 
    double bonds to atoms at the same time, because problems. So,
    all double bonds (or single bonds, or / bonds, or ...) which precede
    a letter shall be replaced with a SPECIAL UNIQUE SEQUENCE.
    For now, that unique sequence is "the same thing, but surrounded by
    tildes".
    """
    output = smiles
    for key, val in FOR_PREPROCESSOR.items():
        output = re.sub(
            key + r'([' + ATOM_MATCHER + r'])',
            val + r'\1',
            output
        )
    return output


##################
##### PARSER #####
##################

PG = ParserGenerator(
    SYMBOLS+LETTERS+['DIGIT', 'TERMINATOR']+BOND_SYMBOLS+BOND_SYMBOLS_TILDE,
    precedence=[],
    cache_id='molparser',
)

# main :: [Molecule].
@PG.production("main : smiles")
def main_production(p):
    "Return the result."
    return p[0]

# smiles ::= terminator | chain terminator
# smiles :: [Molecule].
@PG.production("smiles : chain terminator")
@PG.production("smiles : chain")
def smiles_production(p):
    "Return a list of molecules."
    chain = p[0]
    assert_isinstance(chain, Chain)

    ## Postprocessing!

    ## :: [(Atom, int ring_index, str bond_char)]. Order matters!
    ring_data = chain.ring_data_list

    molecule = chain.molecule
    pending_rings = {}  ## { int ring_index : (atom, str bond_char)}
    for atom, ring_index, bond in ring_data:
        if ring_index in pending_rings.keys():
            other_atom, other_bond = pending_rings[ring_index]
            del pending_rings[ring_index]
            if not (other_bond == DEFAULT_BOND or bond == DEFAULT_BOND):
                assert other_bond == bond, \
                "Unequal ringlinks used: %s, %s" % (other_bond, bond)
            bond_to_use = other_bond if bond == DEFAULT_BOND else bond
            bond_at(bond_to_use, molecule, other_atom, atom, molecule)
        else:
            pending_rings[ring_index] = (atom, bond)

    assert len(pending_rings.keys()) == 0, \
        "Invalid! Unfinished ring bonds! %s" % repr(pending_rings.keys())

    ## TODO: Post-processing of aromaticity

    ## TODO: Post-processing of chirality

    chain.molecule.addHydrogens()
    for mol in chain.dotted_molecules:
        mol.addHydrogens()

    return [chain.molecule] + chain.dotted_molecules

@PG.production("smiles : terminator")
def smiles_empty(p):
    "Return an empty list of molecules."
    return []

# terminator ::= SPACE TAB | LINEFEED | CARRIAGE_RETURN | END_OF_STRING
# terminator :: None.
@PG.production("terminator : TERMINATOR")
def terminator_production(p):
    "Terminators are unimportant."
    return None


##### BONDS #####
# http://www.opensmiles.org/opensmiles.html#bonds

class Chain(BaseBox):
    "Represents a node in the AST for the `chain` variable."
    def __init__(self, mol, first, last):
        """mol :: Molecule.
        first, last :: Atom."""
        self.molecule = mol
        self.first_atom = first
        self.last_atom = last
        self.dotted_molecules = []
        self.ring_data_list = []

    def append_dotteds(self, molecules):
        "molecule :: [Molecule]."
        for molecule in molecules:
            self.dotted_molecules.append(molecule)

    # def append_ring(self, ring_index, bond_char):
    #     """ring_index :: int.
    #     bond_char :: str."""
    #     self.ring_data_list.append((self, ring_index, bond_char))

    def append_rings(self, ring_list):
        """ring_list :: [(Atom, int ring_index, str bond_char)]."""
        self.ring_data_list += ring_list


# chain ::= branched_atom | chain branched_atom |
#           chain bond branched_atom | chain dot branched_atom
# chain :: Chain.
@PG.production("chain : branched_atom")
def chain_production(p):
    "Convert a branched_atom into a chain."
    mol, atom = p[0].as_tuple()
    assert_isinstance(mol, Molecule)
    assert_isinstance(atom, Atom)
    chain = Chain(mol, atom, atom)
    chain.append_dotteds(p[0].dotted_molecules)
    chain.append_rings(p[0].ring_data_list)
    return chain
@PG.production("chain : chain unparenbranch")
def chain_unparenbranch_production(p):
    "Attach a chain to the following branch (the rest of the chain)."
    chain = p[0]
    branch = p[1]
    assert_isinstance(chain, Chain)
    assert_isinstance(branch, Branch)
    ## Now, ligate branch onto chain at the last atom of chain
    bond_to_use = branch.bond_or_dot
    is_dot = (bond_to_use == '.')
    mol1, atom1 = chain.molecule, chain.last_atom
    mol2, atom2 = branch.chain.molecule, branch.chain.first_atom
    assert_isinstance(mol1, Molecule)
    assert_isinstance(atom1, Atom)
    assert_isinstance(mol2, Molecule)
    assert_isinstance(atom2, Atom)
    if is_dot:
        chain.append_dotteds([mol2])
        chain.append_dotteds(branch.chain.dotted_molecules)
        ## Atom, int, str -> int, char
        chain.append_rings(branch.chain.ring_data_list)
        return chain
    else:
        bond_at(bond_to_use, mol1, atom1, atom2, mol2)
        chain.last_atom = branch.chain.last_atom
        chain.append_dotteds(branch.chain.dotted_molecules)
        chain.append_rings(branch.chain.ring_data_list)
        return chain

# bond ::= '-' | '=' | '#' | '$' | ':' | '/' | '\\'
# bond :: str
@PG.production("bond : -")
@PG.production("bond : =")
@PG.production("bond : #")
@PG.production("bond : $")
@PG.production("bond : colon")
@PG.production("bond : /")
@PG.production("bond : \\")
def bond_production(p):
    "Return the string that the bond token has."
    bond = p[0].getstr()
    assert_isinstance(bond, basestring)
    return bond
@PG.production("letterbond : ~-~")
@PG.production("letterbond : ~=~")
@PG.production("letterbond : ~#~")
@PG.production("letterbond : ~$~")
@PG.production("letterbond : ~colon~")
@PG.production("letterbond : ~/~")
@PG.production("letterbond : ~\\~")
def letterbond_production(p):
    "Return the string that the bond token has. But: no tildes!"
    bond = p[0].getstr()
    assert_isinstance(bond, basestring)
    return bond.replace("~", "")


# dot ::= '.'
# dot :: str
@PG.production("dot : .")
def dot_production(p):
    "Return '.'."
    dot = p[0].getstr()
    assert_isinstance(dot, basestring)
    return dot

class BranchedAtom(BaseBox):
    "A class to represent a node in the AST, for the `branched_atom` variable."
    def __init__(self, mol, atom):
        """mol :: Molecule.
        atom :: Atom."""
        self.molecule = mol
        self.atom = atom
        self.dotted_molecules = []
        ## ring_data_list :: [(Atom, int ring_index, str bond_char)].
        self.ring_data_list = []

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
        self.ring_data_list.append((self.atom, ring_index, bond_char))

    def append_rings(self, ring_list):
        """ring_list :: [(Atom, int ring_index, str bond_char)]."""
        self.ring_data_list += ring_list


def add_rings_and_branches(atom, ringlinks, branches):
    "return :: BranchedAtom"
    assert_isinstance(atom, Atom)
    assert_isinstance(ringlinks, list)
    assert_isinstance(branches, list)
    output = BranchedAtom(Molecule(atom), atom)
    for branch in branches:
        output.append_dotteds(branch.chain.dotted_molecules)
        output.append_rings(branch.chain.ring_data_list)
        if branch.bond_or_dot == '.':
            output.append_dotteds([branch.chain.molecule])
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
    for ringlink in ringlinks:
        output.append_ring(ringlink.index, ringlink.bond)
    return output

# branched_atom :: BranchedAtom.
# middle_atom :: Atom, [Ringlink], [Branch].
# ringed_atom :: Atom, [Ringlink], [Branch].
@PG.production("branched_atom : middle_atom")
def branched_to_branched2(p):
    " :: BranchedAtom"
    atom, ringlinks, branches = p[0]
    return add_rings_and_branches(atom, ringlinks, branches)
@PG.production("middle_atom : ringed_atom")
def branched_to_ringed(p):
    " :: Atom, [Ringlink], [Branch]"
    return p[0]
@PG.production("middle_atom : middle_atom branch")
def branched_to_branched(p):
    " :: Atom, [Ringlink], [Branch]"
    atom, ringlinks, branches = p[0]
    return atom, ringlinks, branches+[p[1]]
@PG.production("ringed_atom : ringed_atom ringlink")
def ringed_to_ringed(p):
    " :: Atom, [Ringlink], [Branch]"
    atom, ringlinks, branches = p[0]
    return atom, ringlinks+[p[1]], branches
@PG.production("ringed_atom : atom")
def ringed_to_atom(p):
    " :: Atom, [Ringlink], [Branch]"
    atom = p[0]
    return atom, [], []


class Ringlink(BaseBox):
    "An AST node to store data for the variable `ringlink`."
    def __init__(self, index, bond):
        assert_isinstance(index, int)
        assert_isinstance(bond, basestring)
        self.index = index
        self.bond = bond

class Branch(BaseBox):
    "An AST node to store data for the variable `branch`."
    def __init__(self, bond_or_dot, chain):
        self.bond_or_dot = bond_or_dot
        self.chain = chain

# branch ::= '(' chain ')' | '(' bond chain ')' | '(' dot chain ')'
# branch :: Branch.
# unparenbranch :: Branch.
@PG.production("branch : ( unparenbranch )")
def branch_production_unparen(p):
    "Extract a branch from within parentheses, and return as-is."
    assert p[0].getstr() == '(' and p[2].getstr() == ')', \
        "Not enclosed in parentheses: %s" % str(p)
    return p[1]
@PG.production("unparenbranch : chain")
def branch_production(p):
    "Convert a chain to a branch: add the implicit default bond."
    chain = p[0]
    return Branch(DEFAULT_BOND, chain)
@PG.production("unparenbranch : letterbond chain")
def branch_bond(p):
    "Convert a chain to a branch: add the specified bond."
    bond = p[0]
    chain = p[1]
    assert_isinstance(bond, basestring)
    assert_isinstance(chain, Chain)
    return Branch(bond, chain)
@PG.production("unparenbranch : dot chain")
def branch_dot(p):
    "Convert a chain to a branch: use the dot as a placeholder bond for now."
    dot = p[0]
    chain = p[1]
    assert_isinstance(dot, basestring)
    assert_isinstance(chain, Chain)
    return Branch(dot, chain)

# ringlink ::= bond? DIGIT | bond? '%' DIGIT DIGIT
# ringlink :: Ringlink
@PG.production("ringlink : ringlinkbond")
def ringlink_to_ringlinkbond(p):
    "Allow a ringlink to become a ringlinkbond, which allows specifying a bond."
    return p[0]
@PG.production("ringlink : DIGIT")
def ringlink_prod_1(p):
    "e.g. CC1CCCCC1"
    return ringlink_prod_2([DEFAULT_BOND] + p)
@PG.production("ringlinkbond : bond DIGIT")
def ringlink_prod_2(p):
    "e.g. CCC=1CCC=1"
    bond = p[0]
    index = int(p[1].getstr())
    return Ringlink(index, bond)
@PG.production("ringlink : percentdigit")
def ringlink_prod_3(p):
    "e.g. C%11CCCC%10CC"
    return Ringlink(p[0], DEFAULT_BOND)
@PG.production("ringlinkbond : bond percentdigit")
def ringlink_prod_4(p):
    "e.g. C=%10CCC%10"
    bond = p[0]
    return Ringlink(p[1], bond)

# percentdigit : int
@PG.production("percentdigit : % DIGIT DIGIT")
def percentdigit(p):
    "Convert e.g. '%99' to e.g. int(99)."
    first = p[0]
    assert first.getstr() == '%', "Problem: %s instead of %s?"%(first, '%')
    index_1 = int_from_digit(p[1])
    index_2 = int_from_digit(p[2])
    index = 10*index_1 + index_2
    return index

##### ATOMS #####

# atom ::= bracket_atom | aliphatic_organic | aromatic_organic | '*'
# atom :: Atom.
@PG.production("atom : bracket_atom")
@PG.production("atom : aliphatic_organic")
@PG.production("atom : aromatic_organic")
def atom_production(p):
    "return :: Atom."
    return p[0]
@PG.production("atom : *")
def atom_production_wildcard(p):
    "return :: Atom. * is the wildcard."
    return Atom("*")

##### ORGANIC SUBSET ATOMS #####

# aliphatic_organic ::= 'B'| 'C'| 'N'| 'O'| 'S'| 'P'| 'F'| 'Cl'| 'Br'| 'I'
# aliphatic_organic :: Atom.
ALIPHATIC_ORGANIC = ['B', 'C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I']
def make_aliphatic_organic_production(aliphatic):
    @PG.production("aliphatic_organic : %s" % (' '.join(aliphatic)))
    def _(p):
        "return :: Atom."
        return Atom(''.join(i.getstr() for i in p))
for aliphatic in ALIPHATIC_ORGANIC:
    make_aliphatic_organic_production(aliphatic)

# aromatic_organic ::= 'b' | 'c' | 'n' | 'o' | 's' | 'p'
# aromatic_organic :: Atom.
AROMATIC_ORGANIC = ['b', 'c', 'n', 'o', 's', 'p']
def make_aromatic_organic_production(aromatic):
    @PG.production("aromatic_organic : %s" % (' '.join(aromatic)))
    def _(p):
        "return :: Atom. Also sets is_aromatic to True for future processing."
        assert_isinstance(p[0].getstr(), basestring)
        output = Atom(p[0].getstr().upper())
        output.is_aromatic = True
        return output
for aromatic in AROMATIC_ORGANIC:
    make_aromatic_organic_production(aromatic)

##### BRACKET ATOMS #####

bracket_atom_production_formula = "bracket_atom : [ %s symbol %s %s %s %s ]"
bracket_atom_parts = ["isotope", "chiral", "hcount", "charge", "class"]

def make_bracket_production(five_bools):
    formats = [bracket_atom_parts[i] if five_bools[i] else "" for i in range(5)]
    @PG.production(bracket_atom_production_formula % tuple(formats))
    def _(p):
        "Produce a bracket_atom from some maybe-present parameters."
        six_bools = five_bools[:1] + (True,) + five_bools[1:]
        k = 1 ## start AFTER the [
        args = []
        for j in range(6):
            if six_bools[j]:
                args.append(p[k])
                k += 1
            else:
                args.append(None)
        return bracket_atom_full(*args)
for five_bools in itertools.product((True, False), repeat=5):
    make_bracket_production(five_bools)

def bracket_atom_full(isot, symb, chir, hcou, chge, clss):
    """
    Produce an Atom from various parameters.
    isot :: int.
    symb :: str.
    chir :: Chirality.
    hcou :: int.
    chge :: int.
    clss :: int.
    return :: Atom.
    """
    output = Atom(symb)
    output.isotope = isot
    output.chirality = chir
    output.hcount = 0 if (hcou is None) else hcou
    output.charge = 0 if (chge is None) else chge
    # output.atom_class = clss   ## "Atom class" is throwaway information.

    ## TODO: Check if aromatic. Flag is_aromatic if so.
    return output


# symbol ::= element_symbols | aromatic_symbols | '*'
# symbol :: str.
@PG.production("symbol : *")
def symbol_production_from_wildcard(p):
    "return :: str. * is the wildcard."
    return "*"
@PG.production("symbol : element_symbols")
@PG.production("symbol : aromatic_symbols")
def symbol_production_from_element(p):
    "return :: str."
    return p[0]

# isotope ::= NUMBER
# isotope :: int.
@PG.production("isotope : NUMBER")
def isotope_production(p):
    "return :: int."
    return p[0]

# NUMBER :: int.
@PG.production("NUMBER : NUMBER DIGIT")
def number_recursive_production(p):
    "return :: int."
    return p[0] * 10 + int_from_digit(p[1])
@PG.production("NUMBER : DIGIT")
def number_base_production(p):
    "return :: int."
    num = int_from_digit(p[0])
    return num

# element_symbols ::= [THE ENTIRE PERIODIC TABLE]
# element_symbols :: str.
ELEMENT_SYMBOLS = [
    'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg',
    'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
    'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb',
    'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In',
    'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'Hf', 'Ta', 'W', 'Re', 'Os',
    'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra',
    'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Fl', 'Lv', 'La',
    'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm',
    'Yb', 'Lu', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf',
    'Es', 'Fm', 'Md', 'No', 'Lr']
def make_element_production(element):
    @PG.production("element_symbols : %s" % (' '.join(element)))
    def _(p):
        "return :: str."
        return ''.join(i.getstr() for i in p)
for element in ELEMENT_SYMBOLS:
    make_element_production(element)

# aromatic_symbols ::= 'b' | 'c' | 'n' | 'o' | 'p' | 's' | 'se' | 'as'
# aromatic_symbols :: str.
AROMATIC_SYMBOLS = ['b', 'c', 'n', 'o', 'p', 's', 'se', 'as']
def make_aromatic_production(aromatic):
    @PG.production("aromatic_symbols : %s" % (' '.join(aromatic)))
    def _(p):
        "return :: str."
        return ''.join(i.getstr() for i in p)
for aromatic in AROMATIC_SYMBOLS:
    make_aromatic_production(aromatic)

##### CHIRALITY #####

# chiral ::= '@' | '@@' | '@TH1' | '@TH2' | '@AL1' | '@AL2' | '@SP1' | 
#'@SP2' | '@SP3' | '@TB1' | '@TB2' | '@TB3' | ... | '@TB20' | '@OH1' | 
#'@OH2' | '@OH3' | ... | '@OH30' | '@TB' DIGIT DIGIT | '@OH' DIGIT DIGIT
# chiral :: Chirality
## TODO other productions besides @
@PG.production("chiral : @")
@PG.production("chiral : @ @")
def chiral_production(p):
    "return :: Chirality."
    return Chirality(''.join([i.getstr() for i in p]))

class Chirality(BaseBox):
    """
    Stores chirality information.
    Currently useless.
    """
    def __init__(self, string):
        self.string = string

##### HYDROGENS #####

# hcount ::= 'H' | 'H' DIGIT
# hcount :: int
@PG.production("hcount : H")
def hcount_production_single(p):
    "return :: int"
    return 1
@PG.production("hcount : H DIGIT")
def hcount_production_many(p):
    "return :: int"
    return int_from_digit(p[1])

##### CHARGES #####

# charge ::= '-' | '-' DIGIT? DIGIT | '+' | '+' DIGIT? DIGIT | '--' 
#deprecated | '++' deprecated
# charge :: int
@PG.production("charge : -")
def charge_production_single_minus(p):
    return -1
@PG.production("charge : +")
def charge_production_single_plus(p):
    return +1
@PG.production("charge : - -")
def charge_production_double_minus(p):
    return -2
@PG.production("charge : + +")
def charge_production_double_plus(p):
    return +2
@PG.production("charge : - DIGIT")
def charge_production_minus_many(p):
    return -int_from_digit(p[1])
@PG.production("charge : + DIGIT")
def charge_production_plus_many(p):
    return +int_from_digit(p[1])


##### ATOM CLASS #####

# class ::= ':' NUMBER
# class :: int
@PG.production("class : colon NUMBER")
def class_production(p):
    "return :: int."
    return p[1]

@PG.error
def error_handler(token, expected=None):
    "Handle parser errors."
    raise ValueError(("Ran into a %s (%s) where it wasn't expected."+\
        "At %s. Instead expected: %s.") % (repr(token.name), \
        repr(token.value), dictof(token.source_pos), repr(expected)))


LEXER = LG.build()
PARSER = PG.build()






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
    r"N1CCN(CC1)C(C(F)=C2)=CC(=C2C4=O)N(C3CC3)C=C4(=O)O",
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
    r"OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H](O)[C@@H](O)1", ## glucose
    r"CC[C@H](O1)CC[C@@]12CCCO2", ## another pheromone
    r"CC(C)[C@@]12C[C@@H]1[C@@H](C)C(=O)C2", ## alpha-thujone
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
example_smiles_hard_no_chirality = [
    r"CCC[CH](O)CCC=CC=CC#CC#CC=CCO", ## enanthotoxin
    r"COC(=O)C(C)=CC1C(C)(C)[CH]1C(=O)O[CH]2C(C)=C(C(=O)C2)CC=CC=C",
    r"CC(=O)OCCC(C)=CC[CH](C(C)=C)CCC=C", ## some pheromone 
    r"C[C](C)(O1)C[CH](O)[C]1(O2)[CH](C)[CH]3CC=C4[C]3(C"+\
        "2)C(=O)C[CH]5[CH]4CC[CH](C6)[C]5(C)CC(N7)C6NC(C[C]"+\
        "89(C))C7C[CH]8CC[CH]%10[CH]9C[CH](O)[C]%11(C)C%"+\
        "10=C[CH](O%12)[C]%11(O)[CH](C)[C]%12(O%13)[CH](O)C[C"+\
        "]%13(C)CO", ## Cephalostatin-1 without aromaticity. Note 
                       ## the % before ring closure labels above 9.
]

example_smiles = example_smiles_easy + example_smiles_cistrans +\
    example_smiles_chiral + example_smiles_hard + example_smiles_atoms


class TestToMolecule(unittest.TestCase):
    "Test ALL the features!"

    def tryAssertEqual(self, one, two):
        try:
            self.assertEqual(one, two)
        except AssertionError:
            if DEBUG:
                print "\n%s does not match %s" % (one, two)
                raise

    def assertOne(self, smi):
        "Assert that a single smiles is OK."
        self.tryAssertEqual(
            to_canonical(smi),
            smilesify(moleculify(smi), canonical=True)
        )

    def assertMany(self, smileses):
        "Assert that each in a list of smiles are OK."
        for smi in smileses:
            self.assertOne(smi)

    def assertSame(self, smiles1, smiles2):
        "Assert that two smiles produce the same output."
        self.tryAssertEqual(
            smilesify(moleculify(smiles1), canonical=True),
            smilesify(moleculify(smiles2), canonical=True)
        )
        self.tryAssertEqual(
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
            "C1C(C(C(C)(C)C)C=C(C#C)C2)C2C(C)(C)C=1",
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

    @unittest.skip("Not implemented yet")
    def test_examples_chiral(self):
        "Test all the advanced example smiles."
        smileses = example_smiles_chiral
        self.assertMany(smileses)

    @unittest.skip("Not implemented yet")
    def test_examples_hard(self):
        "Test all the advanced example smiles."
        smileses = example_smiles_hard
        self.assertMany(smileses)

    def test_examples_hard_no_chirality(self):
        "Test all the advanced example smiles."
        smileses = example_smiles_hard_no_chirality
        self.assertMany(smileses)

    def test_examples_atoms(self):
        "Test all the advanced example smiles."
        smileses = example_smiles_atoms
        self.assertMany(smileses)

    def test_wacky_hydrogens(self):
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

    @unittest.skip("Not implemented yet")
    def test_chirality_equivalence(self):
        self.assertSameMany([
            r'N[C@](Br)(O)C',
            r'Br[C@](O)(N)C',
            r'O[C@](Br)(C)N',
            r'Br[C@](C)(O)N',
            r'C[C@](Br)(N)O',
            r'Br[C@](N)(C)O',
            r'C[C@@](Br)(O)N',
            r'Br[C@@](N)(O)C',
            r'[C@@](C)(Br)(O)N',
            r'[C@@](Br)(N)(O)C',
        ])


## TODO temporary
m = moleculify
s = smilesify
