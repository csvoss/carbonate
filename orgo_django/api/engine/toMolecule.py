"""
toMolecule.py

This code contains the method moleculify(smiles), which converts a SMILES string to a Molecule object:
http://en.wikipedia.org/wiki/SMILES

This code is not yet complete.

Only ``moleculify`` is meant to be public-facing. The other methods in this file are private.

This would allow us to take in SMILES strings as inputs, from students or professors, to specify their own molecules to use in problems. Hooray!

Refer to:
    https://github.com/alex/rply
    https://pypi.python.org/pypi/rply/0.5.1
    http://www.opensmiles.org/opensmiles.html
"""

from molecularStructure import *
import string
from rply import ParserGenerator, LexerGenerator
from rply.token import BaseBox


##### HELPER FUNCTIONS #####

def star_it(name, classname):
    """
    Create some productions representing <name>*.
    The productions output a list of <classname> objects.
    name :: str. example: "ringbond"
    classname :: type. example: Ringbond
    starname = name + "star"
    """
    @pg.production("%s : " % starname)
    def star_empty(p):
        return []
    @pg.production("%s : %s %s" % (starname, starname, name))
    def star_full(p):
        xstar = p[0]
        x = p[1]
        assert_isinstance(xstar, list)
        assert_isinstance(x, classname)
        return xstar.append(x)
    return star_empty, star_full

def bond_at(bond, m1, a1, a2, m2):
    BASIC_BONDS = {'-':1, '=':2, '#':3, '$':4}
    if bond == '.':
        raise StandardError("bond_at called with '.' as bond")
    elif bond in BASIC_BONDS.keys():
        ## connect them together
        m1.addMolecule(m2, a2, a1, BASIC_BONDS[bond])
    else:
        raise NotImplementedError(":, /, and \\ are not-yet-implemented bond types.")
        ## TODO.
        ## Idea: Implement / and \\ by adding them as relevant flags,
        ## then postprocessing molecules to convert them into cis/trans centers

def assert_isinstance(i, t):
    assert isinstance(i, t), "In parser: Object %s is of type %s; should be %s" % (repr(i), repr(type(i)), repr(t))


##### LEXER #####

lg = LexerGenerator()
lg.ignore(r"\s+")
#lg.add('SYMBOL', r'[@#$%\*\(\)\[\]=\+\-:/\\\.]') ## SPLIT OUT this one
#lg.add('LETTER', r'[A-IK-PRSTXYZa-ik-pr-vy]') ## SPLIT OUT this one
SYMBOLS = [c for c in "@#$%*()[]=+-:/\\."]
LETTERS = [c for c in "ABCDEFGHIKLMNOPRSTXYZabcdefghiklmnoprstuvy"] # omissions intentional
for sym in SYMBOLS:
    lg.add(sym, sym)
for sym in LETTERS:
    lg.add(sym, sym)

lg.add('DIGIT', r'[0-9]') ## KEEP this one
lg.add('TERMINATOR', r'[ \t\r\n]*') ## KEEP this one




##### PARSER #####
# www.opensmiles.org/opensmiles.html

pg = ParserGenerator(SYMBOLS + LETTERS + ['DIGIT', 'TERMINATOR'], cache_id='molparser')

# main :: [Molecule].
@pg.production("main : smiles")
def main_production(p):
    return p[0]

# smiles ::= terminator | chain terminator
# smiles :: [Molecule].
@pg.production("smiles : chain terminator")
def smiles_production(p):
    chain = p[0]
    assert_isinstance(chain, Chain)
    ## TODO: Post-processing of ringbond list at chain.ring_data_list
    ## TODO: Post-processing of aromaticity
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
        self.ring_data_list = [] ## :: [(Atom, int ring_index, str bond_char)].

    def append_dotteds(self, molecules):
        "molecule :: [Molecule]."
        for molecule in molecules:
            self.dotted_molecules.append(molecule)

    def append_rings(self, ring_index, bond_char):
        """ring_index :: int.
        bond_char :: str."""
        self.ring_data_list.append((self, ring_index, bond_char))

# chain ::= branched_atom | chain branched_atom | chain bond branched_atom | chain dot branched_atom
# chain :: Chain.
@pg.production("chain : branched_atom")
def chain_production(p):
    m, a = p[0].as_tuple()
    assert_isinstance(m, Molecule)
    assert_isinstance(a, Atom)
    return Chain(m, a, a)
@pg.production("chain : chain branched_atom")
def chain_chain(p):
    return chain_bond([p[0], '-', p[1]])
    ##TODO: UNLESS end of chain and atom of branched_atom are aromatic.
    ##Then, we need to add a 1.5 bond there, instead.
    ##And we must be able to tell bond_at that this is the case.
@pg.production("chain : chain bond branched_atom")
def chain_bond(p):
    bond = p[1]
    assert_isinstance(bond, str)
    chain = p[0]
    m1, a1 = chain.molecule, chain.last_atom
    m2, a2 = p[2].as_tuple()
    assert_isinstance(m1, Molecule)
    assert_isinstance(a1, Atom)
    assert_isinstance(m2, Molecule)
    assert_isinstance(a2, Atom)
    bond_at(bond, m1, a1, a2, m2)
    chain.last_atom = a2
    chain.append_dotteds(branched_atom.dotted_molecules)
    chain.append_rings(branched_atom.ring_data_list)
    return chain
@pg.production("chain : chain dot branched_atom")
def chain_dot(p):
    chain = p[0]
    m1, a1 = chain.molecule, chain.last_atom
    dot = p[1]
    m2, a2 = p[2].as_tuple()
    assert dot.getstr() == '.'
    assert_isinstance(m1, Molecule)
    assert_isinstance(a1, Atom)
    assert_isinstance(m2, Molecule)
    assert_isinstance(a2, Atom)
    chain.append_dotteds([m2])
    chain.append_dotteds(branched_atom.dotted_molecules)
    chain.append_rings(branched_atom.ring_data_list)
    return chain


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

    def append_rings(self, ring_index, bond_char):
        """ring_index :: int.
        bond_char :: str."""
        self.ring_data_list.append((self, ring_index, bond_char))


# branched_atom ::= atom ringbond* branch*
# branched_atom :: BranchedAtom
@pg.production("branched_atom : atom ringbondstar branchstar")
def branched_atom(p):
    atom = p[0]
    ringbonds = p[1]
    branches = p[2]
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
        output.append_rings(ringbond.index, ringbond.bond)
    return output

# ringbondstar :: [Ringbond]
ringbondstar_empty, ringbondstar = star_it("ringbond", Ringbond)
# branchstar :: [Branch]
branchstar_empty, branchstar = star_it("branch", Branch)

class Branch(BaseBox):
    def __init__(self, bond_or_dot, chain):
        self.bond_or_dot = bond_or_dot
        self.chain = chain

# branch ::= '(' chain ')' | '(' bond chain ')' | '(' dot chain ')'
# branch :: Branch.
@pg.production("branch : ( chain )")  ## TODO: Make ( and ) be lexable tokens in and of themselves -- as with all symbols, really
def branch_production(p):
    assert p[0]=='(' and p[2]==')'
    chain = p[1]
    return Branch('-', chain)
@pg.production("branch : ( bond chain )")
def branch_bond(p):
    assert p[0]=='(' and p[3]==')'
    bond = p[1]
    chain = p[2]
    assert_isinstance(bond, str)
    assert_isinstance(chain, Chain)
    return Branch(bond, chain)
@pg.production("branch : ( dot chain )")
def branch_dot(p):
    assert p[0]=='(' and p[3]==')'
    dot = p[1]
    chain = p[2]
    assert_isinstance(dot, str)
    assert_isinstance(chain, Chain)
    return Branch(dot, chain)

# bond ::= '-' | '=' | '#' | '$' | ':' | '/' | '\\'
# bond :: str
@pg.production("bond : -")
@pg.production("bond : =")
@pg.production("bond : #")
@pg.production("bond : $")
@pg.production("bond : :")
@pg.production("bond : /")
@pg.production("bond : \\")
def bond_production(p):
    assert_isinstance(p[0], str)
    return p[0]

# dot ::= '.'
# dot :: str
@pg.production("dot : .")
def dot_production(p):
    assert_isinstance(p[0], str)
    return p[0]

class Ringbond(BaseBox):
    def __init__(self, index, bond):
        assert_isinstance(index, int)
        assert_isinstance(bond, str)
        self.index = index
        self.bond = bond

# ringbond ::= bond? DIGIT | bond? '%' DIGIT DIGIT
# ringbond :: Ringbond
@pg.production("ringbond : DIGIT")
def ringbond_prod_1(p):
    ## Again with the TODO about aromatic bonds.
    return ringbond_prod_2(['-'] + p)
@pg.production("ringbond : bond DIGIT")
def ringbond_prod_2(p):
    bond = p[0]
    index = int(p[1])
    return Ringbond(index, bond)
@pg.production("ringbond : % DIGIT DIGIT")
def ringbond_prod_3(p):
    ## A third time with the TODO about aromatic bonds.
    return ringbond_prod_4(['-'] + p)
@pg.production("ringbond : bond % DIGIT DIGIT")
def ringbond_prod_4(p):
    bond = p[0]
    useless_percent = p[1]
    index_1 = int(p[2])
    index_2 = int(p[3])
    index = 10*index_1 + index_2
    return Ringbond(index, bond)

##### ATOMS #####

# atom ::= bracket_atom | aliphatic_organic | aromatic_organic | '*'

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















################ ALL STUFF BELOW IS OLD #############################


@pg.production("expr : expr PLUS expr")
@pg.production("expr : expr MINUS expr")
def expr_op(p):
    lhs = p[0].getint()
    rhs = p[2].getint()
    if p[1].gettokentype() == "PLUS":
        return BoxInt(lhs + rhs)
    elif p[1].gettokentype() == "MINUS":
        return BoxInt(lhs - rhs)
    else:
        raise AssertionError("This is impossible, abort the time machine!")

@pg.production("expr : NUMBER")
def expr_num(p):
    return BoxInt(int(p[0].getstr()))

lexer = lg.build()
parser = pg.build()

class BoxInt(BaseBox):
    def __init__(self, value):
        self.value = value

    def getint(self):
        return self.value


















def moleculify(smiles):
    """
    smiles :: str or [str]. SMILES string(s) e.g. "CC(CN)CCC(O)O"
    return :: Molecule or [Molecule]. Correspondingly.
              If smiles is empty, None.

    Raises a StandardError if the SMILES string contains
    as-yet-unsupported features (like delocalization).
    """

    for i in 'cosn':
        if i in smiles:
            raise StandardError("Unsupported: %s" % i)

    if isinstance(smiles, basestring):
        smiles = str(smiles)
        tokens = lexerSingle(smiles)
        parsed = parserSingle(tokens)
    elif type(smiles) == list:
        tokens = map(lexerSingle, smiles)
        parsed = map(parserSingle, tokens)
    else:
        raise StandardError("Invalid type of input to moleculify: input is %s, type is %s" % (repr(smiles), str(type(smiles))))
    return parsed

    ##return returnExampleMolecule().withHydrogens()

def isAtomString(token):
    """
    token :: str. e.g. Au, [Cu+2], [C@@H]
    return :: bool.
    TEMPORARY -- TODO
    """
    if token in '-=#$/\\':
        return False
    return True

def toAtom(token):
    """
    token :: str. e.g. Au, [Cu+2], [C@@H]
    return :: Atom.
    TEMPORARY -- TODO
    """
    if token[0] != '[':
        return Atom(token)
    else:
        if token[2].lower() in string.lowercase:
            return Atom(token[1:3])
        else:
            return Atom(token[1:2])

def parserSingle(tokens):
    """
    tokens :: [SMILES token (str)].
    return :: Molecule or [Molecule].
              Will only return a list of Molecules if the SMILES is ill-formed
              in such a way as to imply multiple molecules.
    """

    if len(tokens) == 0:
        return None, None

    first = tokens[0]

    if isinstance(first, list):
        ## This... isn't supposed to happen. But it can if the SMILES is weird.
        ## If it does, we treat the initial branch as a *separate molecule*.
        ## e.g. (CCCCC)N(C)Br --> CCCCC and CNBr
        prev_mol, prev_atom = parserWithPassback(first)
        output = [prev_mol]
        molecule = parserSingle(tokens[1:]) ## Recursion!
        if isinstance(molecule, list):
            output += molecule
        elif isinstance(molecule, Molecule):
            output.append(molecule)
        else:
            raise StandardError("parserSingle returned unexpected type: %s, type %s" % (str(molecule), str(type(molecule))))
        return output

    assert isAtomString(first)
    molecule, _ = parserWithPassback(tokens)

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

    ## Comment this to actually test parserSingle
    ## molecule = mol4

    return molecule.withHydrogens()

def parserWithPassback(tokens):
    ## TODO cis trans
    ## TODO chirality
    ## TODO bond types
    ## TODO ring portals
    """
    tokens :: [SMILES token (str)].
    return :: ( Molecule, Atom ).
    """
    first = tokens.pop(0)
    assert isAtomString(first)
    first_atom = toAtom(first)

    molecule = Molecule(first_atom)

    if len(tokens) > 0:
        next = tokens[0]
        while isinstance(next, list):
            next_mol, next_atom = parserWithPassback(next)
            ## attach next_mol to molecule via next_atom connecting to first_atom
            molecule.addMolecule(next_mol, next_atom, first_atom, 1)
            tokens = tokens[1:]
            next = tokens[0]

        if len(tokens) > 0:
            rest_mol, rest_atom = parserWithPassback(tokens)
            molecule.addMolecule(rest_mol, rest_atom, first_atom, 1)

    connecting_atom = first_atom

    return molecule, connecting_atom
    ## TODO: finish this method

def lexerSingle(smiles):
    """
    smiles :: SMILES str.
    return :: [SMILES tokens (str)].
    """
    output, _ = lexerWithPassback(smiles)
    assert _ == None
    return output

def lexerWithPassback(smiles):
    """
    Passback: If you got to this SMILES string by taking off an open-paren,
              this method will tokenize the stuff inside the paren, and return
              both that list of tokens and a 'passback' -- i.e. the remainder
              of the SMILES string.
    smiles :: SMILES str.
    return :: ( [SMILES tokens (str)], remaining SMILES (str) ).
    """
    tokens = []
    while len(smiles) > 0:
        prev_smiles = smiles
        first = smiles[0]
        if first == '[':
            end = smiles.find(']')
            if end == -1:
                raise StandardError("Invalid SMILES: no matching ] in %s" % smiles)
            tokens.append(smiles[0:end+1])
            smiles = smiles[end+1:]
        
        elif first == '(':
            subtree, remaining = lexerWithPassback(smiles[1:])
            tokens.append(subtree)
            smiles = remaining
        
        elif first in '%0123456789':
            if first == '%':
                num, remaining = separateFirstInts(smiles[1:])
            else:
                num = smiles[0]
                remaining = smiles[1:]
            tokens.append(num)
            smiles = remaining
            
        elif first.lower() in string.lowercase:
            ## Check if the next one is a letter too
            try:
                second = smiles[1]
            except IndexError:
                second = '$' ##any non-lowercase char
            if second in string.lowercase:
                tokens.append(first+second)
                smiles = smiles[2:]
            else:
                tokens.append(first)
                smiles = smiles[1:]
        
        elif first in '-=#$/\\':
            tokens.append(first)
            smiles = smiles[1:]

        elif first == ')':
            smiles = smiles[1:]
            return tokens, smiles

        ## no infinite loop!
        assert len(smiles) < len(prev_smiles), "Infinite loop in toMolecule.lexer, unrecognized char: %s" % repr(first)

    assert len(smiles) == 0, "Smiles not fully reduced in toMolecule.lexer"

    return tokens, None


def returnExampleMolecule():
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

    return mol4

def returnExampleSmiles():
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