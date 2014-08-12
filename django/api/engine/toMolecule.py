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
import itertools
import re

from rply import ParserGenerator, LexerGenerator, ParsingError
from rply.token import BaseBox

from molecularStructure import Molecule, Atom, DEBUG
from toCanonical import to_canonical

#### TODO ----
# Tetrahedral Allene-like Systems
# Extended tetrahedral configurations can be specified for conjugated
# allenes with an even number of double bonds. The normal tetrahedral
# rules using '@' and '@@' apply, but the "neighbor" atoms to which the
# chirality refers are at the ends of the allenal system. For example:


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
        return [m for i in smiles for m in moleculify(i)]
    else:
        try:
            if not DEBUG:
                smiles = to_canonical(smiles)
            smiles = preprocess(smiles)
            lexed = LEXER.lex(smiles)
            return PARSER.parse(lexed)
        except ParsingError as e:
            if DEBUG:
                raise StandardError(e.getsourcepos())
            else:
                return Molecule(Atom("*"))


############################
##### HELPER FUNCTIONS #####
############################

BASIC_BONDS = {'-':1, '=':2, '#':3, '$':4}

DEFAULT_BOND = 'default'

## TODO: This is really sketchy and won't work for everything.
## But it works well enough for benzene.
AROMATIC_BOND = 1.5

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
        if DEBUG:
            raise StandardError("bond_at called with '.' as bond")
        else:
            return bond_at(DEFAULT_BOND, mol1, atom1, atom2, mol2)
    elif bond == DEFAULT_BOND:
        ## TODO: support aromaticity!!!
        ## If the two atoms in question are aromatic, this is an aromatic bond
        if atom1.is_aromatic and atom2.is_aromatic:
            mol1.addMolecule(mol2, atom2, atom1, AROMATIC_BOND)
        else:
            mol1.addMolecule(mol2, atom2, atom1, 1)
    elif bond in BASIC_BONDS.keys():
        ## connect them together
        mol1.addMolecule(mol2, atom2, atom1, BASIC_BONDS[bond])
    elif bond == ':':
        ## This is always an aromatic bond.
        mol1.addMolecule(mol2, atom2, atom1, AROMATIC_BOND)
    else:
        ## TODO.
        ## Idea: Implement / and \\ by adding them as relevant flags,
        ## then postprocessing molecules to add cis/trans centers
        if DEBUG:
            raise NotImplementedError("/ and \\ are not yet implemented.")
        else:
            return bond_at(DEFAULT_BOND, mol1, atom1, atom2, mol2)

def dictof(thing):
    "Printable `thing.__dict__` if possible, else just `thing`"
    try:
        return str(thing.__dict__)
    except AttributeError:
        return str(thing)

def debug_decorator(custom_message):
    """
    Apply this decorator to a function to make that function print at you
    whenever it is called.
    """
    def old_debug_decorator(func):
        def new_f(*args, **kwargs):
            "The function, plus a print statement"
            print custom_message
            print "----- Calling %s with %s and %s" % (func, args, kwargs)
            return func(*args, **kwargs)
        return new_f
    return old_debug_decorator

@debug_decorator("Hello!")
def test(x, y):
    print "%s plus %s is %s" % (x, y, str(x+y))

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

SYMBOLS = { ## key (name of token)
    '@': r'@',   ##  : value (regex for lexer to search)
    '%': r'%',
    '*': r'\*',
    '(': r'\(',
    ')': r'\)',
    '[': r'\[',
    ']': r'\]',
    '+': r'\+',
    '.': r'\.',
}
for token_name, regex in SYMBOLS.iteritems():
    LG.add(token_name, regex)

ATOM_MATCHER = r'ABCDEFGHIKLMNOPRSTUVWXYZabcdefghiklmnoprstuvy'
LETTERS = [letter for letter in ATOM_MATCHER]
for token_name in LETTERS:
    LG.add(token_name, token_name)

BOND_SYMBOLS = { ## key (name of token)
    '-': r'\-',       ##  : value (regex for lexer to search)
    '=': r'=',
    '#': r'#',
    '$': r'\$',
    'colon': r':',
    '/': r'/',
    '\\': r'\\',
}
for token_name, regex in BOND_SYMBOLS.iteritems():
    LG.add(token_name, regex)

BOND_SYMBOLS_TILDE = { ## key (name of token)
    '~-~': r'~\-~',         ##  : value (regex for lexer to search)
    '~=~': r'~=~',
    '~#~': r'~#~',
    '~$~': r'~\$~',
    '~colon~': r'~:~',
    '~/~': r'~/~',
    '~\\~': r'~\\~',
}
for token_name, regex in BOND_SYMBOLS_TILDE.iteritems():
    LG.add(token_name, regex)

LG.add('DIGIT', r'[0-9]')
LG.add('TERMINATOR', r'[ \t\r\n]')

FOR_PREPROCESSOR = { ## key (regex for preprocessor to search)
    r'\-': '~-~',    ##  : value (text to replace it with)
    r'=': '~=~',
    r'#': '~#~',
    r'\$': '~$~',
    r':': '~:~',
    r'/': '~/~',
    r'\\': '~\\~',
}


def preprocess(smiles):
    """
    Okay, here's the deal. We can't have double bonds to rings and 
    double bonds to atoms at the same time, because problems. So,
    all double bonds (or single bonds, or / bonds, or ...) which precede
    a letter shall be replaced with a SPECIAL UNIQUE SEQUENCE.
    For now, that unique sequence is "the same thing, but surrounded by
    tildes". See the BOND_SYMBOLS_TILDE and FOR_PREPROCESSOR stuff above.
    """
    output = smiles
    for key, val in FOR_PREPROCESSOR.items():
        output = re.sub(
            key + r'([\[' + ATOM_MATCHER + r'])',
            val + r'\1',
            output
        )

    return output


##################
##### PARSER #####
##################

PG = ParserGenerator(
    SYMBOLS.keys() + LETTERS + ['DIGIT', 'TERMINATOR'] + BOND_SYMBOLS.keys() \
        + BOND_SYMBOLS_TILDE.keys(),
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

            ## Fix the temporary ring indices present in the chirality info
            if atom.chirality is not None:
                ind = atom.chirality.connected_atoms.index(ring_index)
                atom.chirality.connected_atoms[ind] = other_atom
                atom.chirality.update()
            if other_atom.chirality is not None:
                ind = other_atom.chirality.connected_atoms.index(ring_index)
                other_atom.chirality.connected_atoms[ind] = atom
                other_atom.chirality.update()

        else:
            pending_rings[ring_index] = (atom, bond)

    assert len(pending_rings.keys()) == 0, \
        "Invalid! Unfinished ring bonds! %s" % repr(pending_rings.keys())

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

    ## Tetrahedral Chirality
    ## case: last_atom in chain is chiral
    ##       and it needs to know about first_atom in unparenbranch
    ## TODO
    if chain.last_atom.chirality is not None:
        chain.last_atom.chirality.append(branch.chain.first_atom)

    ## case: first_atom in unparenbranch is chiral
    ##       and it needs to know about last_atom in chain
    ## TODO
    if branch.chain.first_atom.chirality is not None:
        branch.chain.first_atom.chirality.prepend(chain.last_atom)

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
    for ringlink in ringlinks:
        output.append_ring(ringlink.index, ringlink.bond)
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
    return output

# branched_atom :: BranchedAtom.
# middle_atom :: Atom, [Ringlink], [Branch].
# ringed_atom :: Atom, [Ringlink].
@PG.production("branched_atom : middle_atom")
def branched_to_branched2(p):
    " :: BranchedAtom"
    atom, ringlinks, branches = p[0]
    return add_rings_and_branches(atom, ringlinks, branches)
@PG.production("middle_atom : ringed_atom")
def branched_to_ringed(p):
    " :: Atom, [Ringlink], [Branch]"
    atoms, ringlinks = p[0]
    branches = []
    return atoms, ringlinks, branches
@PG.production("middle_atom : middle_atom branch")
def branched_to_branched(p):
    " :: Atom, [Ringlink], [Branch]"
    atom, ringlinks, branches = p[0]
    branch = p[1]
    assert_isinstance(atom, Atom)
    assert_isinstance(ringlinks, list)
    assert_isinstance(branches, list)
    assert_isinstance(branch, Branch)

    ## Tetrahedral Chirality
    ## case: middle_atom is chiral
    ##       and it needs to know about first_atom in branch
    if atom.chirality is not None:
        atom.chirality.append(branch.chain.first_atom)

    ## case: first_atom in branch is chiral
    ##       and it needs to know about middle_atom
    if branch.chain.first_atom.chirality is not None:
        branch.chain.first_atom.chirality.prepend(atom)

    return atom, ringlinks, branches+[branch]
@PG.production("ringed_atom : ringed_atom ringlink")
def ringed_to_ringed(p):
    " :: Atom, [Ringlink]"
    atom, ringlinks = p[0]
    ringlink = p[1]
    assert_isinstance(atom, Atom)
    assert_isinstance(ringlinks, list)
    assert_isinstance(ringlink, Ringlink)

    ## Tetrahedral Chirality
    ## case: ringed_atom is chiral
    ##       and it needs to know about the ringlink
    if atom.chirality is not None:
        atom.chirality.append(ringlink.index)

    return atom, ringlinks+[ringlink]
@PG.production("ringed_atom : atom")
def ringed_to_atom(p):
    " :: Atom, [Ringlink]"
    atom = p[0]
    assert_isinstance(atom, Atom)
    return atom, []


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
    "Create several productions, one for each aliphatic-organic."
    @PG.production("aliphatic_organic : %s" % (' '.join(aliphatic)))
    def _(p):
        "return :: Atom."
        return Atom(''.join(i.getstr() for i in p))
for a in ALIPHATIC_ORGANIC:
    make_aliphatic_organic_production(a)

# aromatic_organic ::= 'b' | 'c' | 'n' | 'o' | 's' | 'p'
# aromatic_organic :: Atom.
AROMATIC_ORGANIC = ['b', 'c', 'n', 'o', 's', 'p']
def make_aromatic_organic_production(aromatic):
    "Create several production rules, one for each aromatic-organic."
    @PG.production("aromatic_organic : %s" % (' '.join(aromatic)))
    def _(p):
        "return :: Atom. Also sets is_aromatic to True for future processing."
        assert_isinstance(p[0].getstr(), basestring)
        output = Atom(p[0].getstr().upper())
        output.is_aromatic = True
        return output
for a in AROMATIC_ORGANIC:
    make_aromatic_organic_production(a)

##### BRACKET ATOMS #####

bracket_atom_production_formula = "bracket_atom : [ %s symbol %s %s %s %s ]"
bracket_atom_parts = ["isotope", "chiral", "hcount", "charge", "class"]

def make_bracket_production(five_bools):
    """
    Create several production rules, one for each combination of bracket atom
    parts. Each one may be present or absent. five_bools is a tuple of booleans
    describing the particular conditions to make a production for.
    """
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
for bools in itertools.product((True, False), repeat=5):
    make_bracket_production(bools)

@PG.production("bracket_atom : D")
@PG.production("bracket_atom : T")
def bracket_atom_hydrogen_isotopes(p):
    "Deuterium and tritium."
    if p[0].getstr() == 'D':
        isotope = 2
    elif p[0].getstr() == 'T':
        isotope = 3
    else:
        if DEBUG:
            raise StandardError("Neither D nor T. How did you get here??")
        else:
            isotope = 1
    return bracket_atom_full(isotope, ('H', False), None, None, None, None)

def bracket_atom_full(isot, symb, chir, hcou, chge, clss):
    """
    Produce an Atom from various parameters.
    isot :: int or None.
    symb :: str symbol, bool is_aromatic.
    chir :: Chirality or None.
    hcou :: int or None.
    chge :: int or None.
    clss :: int or None.
    return :: Atom.
    """
    output = Atom(symb[0])
    
    output.isotope = isot
    output.hcount = 0 if (hcou is None) else hcou
    output.charge = 0 if (chge is None) else chge
    output.is_aromatic = symb[1]
    
    output.chirality = chir
    if chir is not None:
        output.chirality.central_atom = output
        if hcou == 0 or hcou is None:
            pass # intentional
        elif hcou == 1:
            output.chirality.append_hydrogen()
        else:
            output.chirality = None
            if DEBUG:
                raise StandardError("No more than 1 H for a single chiral atom")



    # output.atom_class = clss   ## "Atom class" is throwaway information.
    return output


# symbol ::= element_symbols | aromatic_symbols | '*'
# symbol :: str symbol, bool is_aromatic.
@PG.production("symbol : *")
def symbol_production_from_wildcard(p):
    "return :: str, bool. * is the wildcard."
    return "*", False
@PG.production("symbol : element_symbols")
def symbol_production_from_element(p):
    "return :: str, bool. (False for not aromatic)"
    return p[0], False
@PG.production("symbol : aromatic_symbols")
def symbol_production_from_aromatic(p):
    "return :: str, bool. (True for aromaticity)"
    return p[0], True

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
def make_element_production(element_symbol):
    "Create several productions, one for each element."
    @PG.production("element_symbols : %s" % (' '.join(element_symbol)))
    def _(p):
        "return :: str."
        return ''.join(i.getstr() for i in p)
for element in ELEMENT_SYMBOLS:
    make_element_production(element)

# aromatic_symbols ::= 'b' | 'c' | 'n' | 'o' | 'p' | 's' | 'se' | 'as'
# aromatic_symbols :: str.
AROMATIC_SYMBOLS = ['b', 'c', 'n', 'o', 'p', 's', 'se', 'as']
def make_aromatic_production(aromatic):
    "Create several productions, one for each aromatic symbol."
    @PG.production("aromatic_symbols : %s" % (' '.join(aromatic)))
    def _(p):
        "return :: str."
        return ''.join(i.getstr() for i in p)
for a in AROMATIC_SYMBOLS:
    make_aromatic_production(a)

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
    """
    def __init__(self, string, central_atom=None):
        self.string = string
        self.central_atom = central_atom
        self.connected_atoms = []
        self.hydrogen = False

    def append_hydrogen(self):
        "For [C@H] or [C@@H] situations."
        self.connected_atoms.append(None) ## placeholder
        self.hydrogen = True
        self.update()

    def append(self, atom):
        """
        atom :: Atom or int.
        Add atom as a connected element lexically after this chiral-atom.
        """
        # assert_isinstance(atom, Atom)
        self.connected_atoms.append(atom)
        self.update()

    def prepend(self, atom):
        """
        atom :: Atom or int.
        Add atom as a connected element lexically before this chiral-atom.
        """
        # assert_isinstance(atom, Atom)
        self.connected_atoms.insert(0, atom)
        self.update()

    def update(self):
        "Uses known information to update the chirality of central_atom."
        if len(self.connected_atoms) >= 4:
            ref = self.connected_atoms[0]
            cwlist = self.connected_atoms[1:4]
            if self.string == '@': ## counterclockwise
                ccwlist = reversed(cwlist)
                self.central_atom.newChiralCenter(ref, ccwlist)
            elif self.string == '@@': ## clockwise
                self.central_atom.newChiralCenter(ref, cwlist)
            else:
                if DEBUG:
                    raise NotImplementedError("TODO: other forms of chirality")
                else:
                    pass

##### HYDROGENS #####

# hcount ::= 'H' | 'H' DIGIT
# hcount :: int
@PG.production("hcount : H")
def hcount_production_single(p):
    "return :: int"
    return 1
@PG.production("hcount : H NUMBER") ##TODO: spec says DIGIT, but Babel does NUM
def hcount_production_many(p):
    "return :: int"
    # return int_from_digit(p[1])
    return p[1]

##### CHARGES #####

# charge ::= '-' | '-' DIGIT? DIGIT | '+' | '+' DIGIT? DIGIT | '--' 
#deprecated | '++' deprecated
# charge :: int
@PG.production("charge : -")
def charge_production_single_minus(p):
    "-"
    return -1
@PG.production("charge : +")
def charge_production_single_plus(p):
    "+"
    return +1
@PG.production("charge : - -")
def charge_production_double_minus(p):
    "--"
    return -2
@PG.production("charge : + +")
def charge_production_double_plus(p):
    "++"
    return +2
@PG.production("charge : - DIGIT")
def charge_production_minus_many(p):
    "-D"
    return -int_from_digit(p[1])
@PG.production("charge : + DIGIT")
def charge_production_plus_many(p):
    "+D"
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
    if DEBUG:
        raise ValueError(("Ran into a %s (%s) where it wasn't expected."+\
            "At %s. Instead expected: %s.") % (repr(token.name), \
            repr(token.value), dictof(token.source_pos), repr(expected)))
    else:
        print "Warning: parser error"


LEXER = LG.build()
PARSER = PG.build()
