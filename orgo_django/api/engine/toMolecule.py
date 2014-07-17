"""
toMolecule.py

This code contains the method moleculify(smiles), which converts a SMILES string to a Molecule object:
http://en.wikipedia.org/wiki/SMILES

This code is not yet complete. #TODO #YOLO

Only ``moleculify`` is meant to be public-facing. The other methods in this file are private.

This would allow us to take in SMILES strings as inputs, from students or professors, to specify their own molecules to use in problems. Hooray!

"""

from molecularStructure import *
import string

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
        tokens = lexerMany(smiles)
        parsed = parserMany(tokens)
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

def parserMany(tokensMany):
    """
    tokensMany :: [[SMILES token (str)]].
    return :: [Molecule].
    """
    return map(parserSingle, tokensMany)

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

def lexerMany(smiles_list):
    """
    smiles_list :: [SMILES str].
    return :: [[SMILES tokens (str)]].
    """
    ## Many molecules: return a list
    return map(lexerSingle, smiles_list.split('.'))

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