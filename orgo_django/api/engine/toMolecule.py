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
import re
import unittest

from rply import ParserGenerator, LexerGenerator, ParsingError
from rply.token import BaseBox

from molecularStructure import Molecule, Atom
from toSmiles import smilesify, to_canonical



## TODO Remove old, commented-out code
## TODO Refactor hastily-renamed productions, etc.

DEBUG = True

##### HELPER FUNCTIONS #####

BASIC_BONDS = {'-':1, '=':2, '#':3, '$':4}

DEFAULT_BOND = 'default'

# def plus_it(name, classname):
#     """
#     Create some productions representing <name>+.
#     The productions output a list of <classname> objects.
#     name :: str. example: "ringbond"
#     classname :: type. example: Ringbond
#     plusname = name + "plus"
#     """
#     plusname = name + "plus"
#     @PG.production("%s : %s" % (plusname, name))
#     def plus_empty(p):
#         return []
#     @PG.production("%s : %s %s" % (plusname, plusname, name))
#     def plus_full(p):
#         xplus = p[0]
#         x = p[1]
#         assert_isinstance(xplus, list)
#         assert_isinstance(x, classname)
#         return xplus + [x]
#     return plus_empty, plus_full

def bond_at(bond, mol1, atom1, atom2, mol2):
    """
    Create a bond between molecules mol1 and mol2, at atoms atom1 and atom2,
    using bond.

    REQUIRES that atom1 precede atom2 in the SMILES string.

    THIS METHOD IS INCOMPLETE. #TODO
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


def assert_isinstance(instance, some_type):
    """
    Raise an error if type(instance) is not some_type.
    """
    assert isinstance(instance, some_type), \
        "In parser: Object %s is of type %s; should be %s" % \
        (repr(instance), repr(type(instance)), repr(some_type))

##### LEXER #####

## Create a Lexer, which splits up a string into tokens.
## Tokens may be, for example, SYMBOLs or LETTERs or TERMINATORs or DIGITs
## Whitespace will be ignored.

LG = LexerGenerator()
LG.ignore(r"\s+")

#LG.add('SYMBOL', r'[@#$%\*\(\)\[\]=\+\-:/\\\.]') ## SPLIT OUT this one
SYMBOLS_DICT = {
    '@': '@', '%': '%', '*': r'\*', '(': r'\(', ')': r'\)',
    '[': r'\[', ']': r'\]', '+': r'\+', '.': r'\.'
}
SYMBOLS = SYMBOLS_DICT.keys()
for sym in SYMBOLS:
    LG.add(sym, SYMBOLS_DICT[sym])

#LG.add('LETTER', r'[A-IK-PRSTXYZa-ik-pr-vy]') ## SPLIT OUT this one
LETTERS = ["C", "N", "O", "F"] ## TODO this is for debugging
#LETTERS = [c for c in "ABCDEFGHIKLMNOPRSTXYZabcdefghiklmnoprstuvy"]
# (omissions intentional)
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
ATOM_MATCHER = r'ABCDEFGHIKLMNOPRSTXYZabcdefghiklmnoprstuvy'

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


##### PARSER #####
# www.opensmiles.org/opensmiles.html

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
    ## TODO: Post-processing of ringbond list at chain.ring_data_list

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
                "Unequal ringbonds used: %s, %s" % (other_bond, bond)
            bond_to_use = other_bond if bond == DEFAULT_BOND else bond
            bond_at(bond_to_use, molecule, other_atom, atom, molecule)
        else:
            pending_rings[ring_index] = (atom, bond)

    assert len(pending_rings.keys()) == 0, \
        "Invalid! Unfinished ring bonds! %s" % repr(pending_rings.keys())

    ## TODO: Post-processing of aromaticity
    # print smilesify(chain.dotted_molecules)
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
def chain_unparenbranch_testing_temporary(p): #TODO: Does it work?
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

# @PG.production("chain : chain branched_atom")
# def chain_chain(p):
#     return chain_bond([p[0], DEFAULT_BOND, p[1]])
# @PG.production("chain : chain bond branched_atom")
# def chain_bond(p):
#     bond = p[1]
#     assert_isinstance(bond, basestring)
#     chain = p[0]
#     mol1, atom1 = chain.molecule, chain.last_atom
#     mol2, atom2 = p[2].as_tuple()
#     assert_isinstance(mol1, Molecule)
#     assert_isinstance(atom1, Atom)
#     assert_isinstance(mol2, Molecule)
#     assert_isinstance(atom2, Atom)
#     bond_at(bond, mol1, atom1, atom2, mol2)
#     chain.last_atom = atom2
#     chain.append_dotteds(p[2].dotted_molecules)
#     chain.append_rings(p[2].ring_data_list)
#     return chain
# @PG.production("chain : chain dot branched_atom")
# def chain_dot(p):
#     chain = p[0]
#     mol1, atom1 = chain.molecule, chain.last_atom
#     dot = p[1]
#     mol2, atom2 = p[2].as_tuple()
#     assert dot == '.'
#     assert_isinstance(mol1, Molecule)
#     assert_isinstance(atom1, Atom)
#     assert_isinstance(mol2, Molecule)
#     assert_isinstance(atom2, Atom)
#     chain.append_dotteds([mol2])
#     chain.append_dotteds(p[2].dotted_molecules)
#     ## Atom, int, str -> int, char
#     chain.append_rings(p[2].ring_data_list)
#     return chain


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
# @PG.production("bond : BOND")
# def bond_production(p):
#     bond = p[0].getstr()
#     assert_isinstance(bond, basestring)
#     return bond

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


def add_rings_and_branches(atom, ringbonds, branches):
    "return :: BranchedAtom"
    assert_isinstance(atom, Atom)
    assert_isinstance(ringbonds, list)
    assert_isinstance(branches, list)
    output = BranchedAtom(Molecule(atom), atom)
    for branch in branches:
        output.append_dotteds(branch.chain.dotted_molecules)
        output.append_rings(branch.chain.ring_data_list)
        if branch.bond_or_dot == '.':
            pass  ## TODO: Is this meant to be blank?
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
@PG.production("branched_atom : branched_atom2")
def branched_to_branched2(p):
    " :: BranchedAtom"
    atom, ringbonds, branches = p[0]
    return add_rings_and_branches(atom, ringbonds, branches)
## TODO I think this rule is unnecessary.
# @PG.production("branched_atom : branched_atom2 unparenbranch") 
# def branched_to_branched2_with_unparenbranch(p):
#     " :: BranchedAtom"
#     atom, ringbonds, branches = p[0]
#     unparenbranch = p[1]
#     branches = branches + [unparenbranch]
#     return add_rings_and_branches(atom, ringbonds, branches)
@PG.production("branched_atom2 : ringed_atom")
def branched_to_ringed(p):
    " :: Atom, [Ringbond], [Branch]"
    return p[0]
# @PG.production("branched_atom2 : atom")
# def branched_to_plain_atom_maybe_this_will_work(p):
#     " :: Atom, [Ringbond], [Branch]"
#     atom = p[0]
#     return atom, [], []
@PG.production("branched_atom2 : branched_atom2 branch")
def branched_to_branched(p):
    " :: Atom, [Ringbond], [Branch]"
    atom, ringbonds, branches = p[0]
    return atom, ringbonds, branches+[p[1]]
@PG.production("ringed_atom : ringed_atom ringbond")
def ringed_to_ringed(p):
    " :: Atom, [Ringbond], [Branch]"
    atom, ringbonds, branches = p[0]
    return atom, ringbonds+[p[1]], branches
@PG.production("ringed_atom : atom")
def ringed_to_atom(p):
    " :: Atom, [Ringbond], [Branch]"
    atom = p[0]
    return atom, [], []


class Ringbond(BaseBox):
    "An AST node to store data for the variable `ringbond`."
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
def unbranch_production(p):
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

# ringbond ::= bond? DIGIT | bond? '%' DIGIT DIGIT
# ringbond :: Ringbond
@PG.production("ringbond : ringbond2")  ## TODO aaaughh
def ringbond_to_ringbond2(p):
    "Allow a ringbond to become a ringbond2, which allows specifying the bond."
    return p[0]
@PG.production("ringbond : DIGIT")
def ringbond_prod_1(p):
    "e.g. CC1CCCCC1"
    return ringbond_prod_2([DEFAULT_BOND] + p)
@PG.production("ringbond2 : bond DIGIT")
def ringbond_prod_2(p):
    "e.g. CCC=1CCC=1"
    bond = p[0]
    index = int(p[1].getstr())
    return Ringbond(index, bond)
@PG.production("ringbond : percentdigit")
def ringbond_prod_3(p):
    "e.g. C%11CCCC%10CC"
    return Ringbond(p[0], DEFAULT_BOND)
# @PG.production("ringbond2 : bond percentdigit")  ## TODO uncomment
# def ringbond_prod_4(p):
#     bond = p[0]
#     return Ringbond(p[1], bond)

# percentdigit : int
@PG.production("percentdigit : % DIGIT DIGIT")
def percentdigit(p):
    "Convert e.g. '%99' to e.g. int(99)."
    first = p[0]
    assert first.getstr() == '%', "Problem: %s instead of %s?"%(first, '%')
    index_1 = int(p[1].getstr())
    index_2 = int(p[2].getstr())
    index = 10*index_1 + index_2
    return index

##### ATOMS #####

# atom ::= bracket_atom | aliphatic_organic | aromatic_organic | '*'
@PG.production("atom : C")
@PG.production("atom : N")
@PG.production("atom : F")
@PG.production("atom : O") ## TODO: atoms better
def atom_production(p):
    'Return an Atom object!'
    ## TODO -- Temporary
    return Atom(p[0].getstr())

##### ORGANIC SUBSET ATOMS #####

# aliphatic_organic ::= 'B'| 'C'| 'N'| 'O'| 'S'| 'P'| 'F'| 'Cl'| 'Br'| 'I'

# aromatic_organic ::= 'b' | 'c' | 'n' | 'o' | 's' | 'p'

##### BRACKET ATOMS #####

# bracket_atom ::= '[' isotope? symbol chiral? hcount? charge? class? ']'

# symbol ::= element_symbols | aromatic_symbols | '*'

# isotope ::= NUMBER

# element_symbols ::= 'H' | 'He' | 'Li' | 'Be' | 'B' | 'C' | 'N' | 'O' |
# 'F' | 'Ne' | 'Na' | 'Mg' | 'Al' | 'Si' | 'P' | 'S' | 'Cl' | 'Ar' | 'K' 
#| 'Ca' | 'Sc' | 'Ti' | 'V' | 'Cr' | 'Mn' | 'Fe' | 'Co' | 'Ni' | 'Cu' | 
#'Zn' | 'Ga' | 'Ge' | 'As' | 'Se' | 'Br' | 'Kr' | 'Rb' | 'Sr' | 'Y' | 
#'Zr' | 'Nb' | 'Mo' | 'Tc' | 'Ru' | 'Rh' | 'Pd' | 'Ag' | 'Cd' | 'In' | 
#'Sn' | 'Sb' | 'Te' | 'I' | 'Xe' | 'Cs' | 'Ba' | 'Hf' | 'Ta' | 'W' | 'Re' 
#| 'Os' | 'Ir' | 'Pt' | 'Au' | 'Hg' | 'Tl' | 'Pb' | 'Bi' | 'Po' | 'At' | 
#'Rn' | 'Fr' | 'Ra' | 'Rf' | 'Db' | 'Sg' | 'Bh' | 'Hs' | 'Mt' | 'Ds' | 
#'Rg' | 'Cn' | 'Fl' | 'Lv' | 'La' | 'Ce' | 'Pr' | 'Nd' | 'Pm' | 'Sm' | 
#'Eu' | 'Gd' | 'Tb' | 'Dy' | 'Ho' | 'Er' | 'Tm' | 'Yb' | 'Lu' | 'Ac' | 
#'Th' | 'Pa' | 'U' | 'Np' | 'Pu' | 'Am' | 'Cm' | 'Bk' | 'Cf' | 'Es' | 
#'Fm' | 'Md' | 'No' | 'Lr'

# aromatic_symbols ::= 'b' | 'c' | 'n' | 'o' | 'p' | 's' | 'se' | 'as'

##### CHIRALITY #####

# chiral ::= '@' | '@@' | '@TH1' | '@TH2' | '@AL1' | '@AL2' | '@SP1' | 
#'@SP2' | '@SP3' | '@TB1' | '@TB2' | '@TB3' | ... | '@TB20' | '@OH1' | 
#'@OH2' | '@OH3' | ... | '@OH30' | '@TB' DIGIT DIGIT | '@OH' DIGIT DIGIT

##### HYDROGENS #####

# hcount ::= 'H' | 'H' DIGIT

##### CHARGES #####

# charge ::= '-' | '-' DIGIT? DIGIT | '+' | '+' DIGIT? DIGIT | '--' 
#deprecated | '++' deprecated

##### ATOM CLASS #####

# class ::= ':' NUMBER



def dictof(thing):
    "Printable `thing.__dict__` if possible, else just `thing`"
    try:
        return str(thing.__dict__)
    except AttributeError:
        return str(thing)

@PG.error
def error_handler(token, expected=None):
    "Handle parser errors."
    raise ValueError(("Ran into a %s (%s) where it wasn't expected."+\
        "At %s. Instead expected: %s.") % (repr(token.name), \
        repr(token.value), dictof(token.source_pos), repr(expected)))



LEXER = LG.build()
PARSER = PG.build()



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
]
example_smiles_atoms = [
    r"[2H]C(Cl)(Cl)Cl",  ## deuterochloroform (hydrogen-2)
    r"[Cu+2].[O-]S(=O)(=O)[O-]" ## copper(ii) sulfate
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

example_smiles = example_smiles_easy + example_smiles_cistrans +\
    example_smiles_chiral + example_smiles_hard + example_smiles_atoms



##### UNIT TESTS #####

class TestToMolecule(unittest.TestCase):
    "Test ALL the features!"

    def assertOkay(self, smi):
        "Assert that a single smiles is OK."
        try:
            self.assertEqual(set([to_canonical(smi)]),
                             set(smilesify(moleculify(smi))))
        except AssertionError:
            if DEBUG:
                print "\nFailed on %s with bad output %s" % \
                    (smi, smilesify(moleculify(smi)))
                raise

    def assertMany(self, smileses):
        "Assert that each in a list of smiles are OK."
        for smi in smileses:
            self.assertOkay(smi)

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
            self.assertOkay(smi)

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

    def test_rings(self):
        "Test that basic usage of ringbonds is functional."
        smileses = [
            r"C1CCCC1",  #TODO uncomment
            r"CCC2CC(CC(CC)(C2)CC)CC",
            r"CCC(CC)(CC1)CCCC1",
            r"C1(C1)",
        ]
        self.assertMany(smileses)

    def test_rings_many(self):
        "Test that multiple ringbonds can be used at the same time."
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
        "Test that percent-escaping ringbonds works."
        smileses = [
            r"CCC%10CCCC%10",
            r"CCC%99CCC%22CCC%22CCCC%99",
            r"CCC1CC%12CCC1CC%12CCC",
        ]
        self.assertMany(smileses)

    def test_rings_with_bonds(self):
        "Test that prefacing ringbonds with bonds works."
        smileses = [
            r"CCC=1CCC1",
            r"CCC=1CCCC=1",
            r"C1CCCC#2CCC1C2C",
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
            first = set([to_canonical(i) for i in split])
            second = set(smilesify(moleculify(dotted)))
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

    @unittest.skip("Not implemented yet")
    def test_examples_atoms(self):
        "Test all the advanced example smiles."
        smileses = example_smiles_atoms
        self.assertMany(smileses)


## TODO temporary
m = moleculify
s = smilesify
