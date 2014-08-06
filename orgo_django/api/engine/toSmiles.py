"""
toSmiles.py

This code contains the function smilesify(molecule), which converts a Molecule
object to a canonical SMILES string:
http://en.wikipedia.org/wiki/SMILES
Once a Molecule has been converted to a SMILES string, we can render it as an
SVG using renderSVG.py.

PROCEED CAUTIOUSLY -- THIS CODE IS FINICKY

Public-facing methods:
    `smilesify`

Does not yet support:
- charge :( :( 
- radicals :( :( 
- isotopes
"""
from molecularStructure import Molecule, Atom
import copy
from toCanonical import to_canonical


def smilesify(molecule, canonical=True):
    """
    molecule :: Molecule or [Molecule]
    (in the latter case, smilesify(mol) is applied for mol in the list)

    Converts Molecule objects to SMILES strings. Traverses the molecule to
    detect and mark rings or cycles. Uses a convoluted system of flags to 
    do this. Then, passes it on to _subsmiles, which operates on a molecule
    with rings already flagged, and performs tree traversal.

    return :: str.
    """
    if isinstance(molecule, list):
        try:
            _assertMolecule(molecule[0])
        except IndexError:
            return ""
        output = '.'.join(
            [smilesify(molec, canonical=canonical) for molec in molecule]
        )
        if canonical:
            return to_canonical(output)
        else:
            return output

    if len(molecule.atoms) == 0:
        return ""

    _initialize_nonHNeighbors(molecule)
    _flag_rings(molecule)
    
    output = _get_generated_smiles(molecule)
    _resetflags(molecule)

    if canonical:
        return to_canonical(output)
    else:
        return output

BOND_SYMBOLS = {0: '0', 1: '-', 2: '=', 3: '#', 4: '$', 1.5: ':'}

def _subsmiles(molecule, start_atom, parent_atom):
    """
    Precondition: molecule has been flagged for ring positioning (some rflag 
        values on atoms might != 0). This is done by smilesify().
    Creates and returns a SMILES string for unflagged (!atom.flag==2) atoms
        within a molecule, starting with the given atom.

    molecule :: Molecule.
    start_atom :: Atom.
    parent_atom :: Atom or 0.

    Traverses the molecule from the given starting atom, returning the SMILES
     representation.
    Can be called recursively -- hence start_atom and parent_atom. This is a
     tree traversal, after all!

    return :: a SMILES substring
    """

    assert isinstance(molecule, Molecule), \
        "molecule invalid: %s" % str(molecule)
    assert isinstance(start_atom, Atom), \
        "start_atom invalid: %s" % str(start_atom)
    try:
        assert isinstance(parent_atom, Atom), "parent_atom invalid"
    except AssertionError:
        assert parent_atom == 0, \
            "parent_atom invalid: %s" % str(parent_atom)
    
    #Flag the current atom.
    start_atom.flag = 2
    
    output = str(start_atom)

    if start_atom.is_cistrans:
        return _get_subsmiles_ct(output, molecule, start_atom, parent_atom)
    

    ###Put a ring marker on the atom, if its ring partner is not flagged yet.
    ##if (start_atom.rflag != 0) and (start_atom.rAtom.flag != 2):
    ##    output += str(start_atom.rflag)

    #Check if the atom is a chiral center. If so:
    if start_atom.is_chiral:
        hasP = (parent_atom != 0)
        hasH = (True in \
                    [a.element.lower() == "h" \
                     for a in list(start_atom.neighbors)]) or \
               (None in \
                    [start_atom.chiralA, start_atom.chiralB,
                     start_atom.chiralC, start_atom.chiralD])
        #If the atom has a hydrogen:
        #Add [ and @@H] to the current output. (e.g. [C@@H]
        #If the atom does not have a hydrogen:
        #Add [ and @@] to the current output. (e.g. [C@@]
        #Then, add the remaining neighbors in the proper order.
        #Implement this by smart rearrangement of to_add.
        #Be mindful of: whether or not there is a parent atom; "" a hydrogen.
        if hasH:
            # output = "[" + output + "@@H]"
            if hasP:
                #to_add should have two elements
                l = start_atom.chiralCWlist(parent_atom) #list of three atoms
                if None in l:
                    x = l.index(None) #index of hydrogen atom in list
                else:
                    x = 0
                    while x < len(l):
                        if l[x].element == 'H':
                            break
                        x += 1
                to_add = [l[(x+1) %3], l[(x+2) %3]] #correct permutation
                if None in to_add:
                    raise StandardError("%s is chiral, but has two hydrogens." \
                        % start_atom.element)
            else:
                #to_add should have three elements
                to_add = start_atom.chiralCWlist(None) #list of three atoms
                if None in to_add:
                    raise StandardError("%s is chiral, but has two hydrogens." \
                        % start_atom.element)
        else:
            # output = "[" + output + "@@]"
            if hasP:
                #to_add should have three elements
                to_add = start_atom.chiralCWlist(parent_atom)
            else:
                #to_add should have four elements
                arbitraryRef = list(start_atom.neighbors)[0]
                l = start_atom.chiralCWlist(arbitraryRef)
                to_add = [arbitraryRef] + l

    
    # Prepare to add new groups for all neighbor atoms which are not the parent
    # atom and not the rAtom.
    else:
        to_add = [atom for atom in list(start_atom.non_h_neighbors) if not \
            (atom == parent_atom or atom == None)]

    # if to_add == None:
    #     to_add = [atom for atom in list(start_atom.non_h_neighbors) if not \
    #           (atom==parent_atom or atom==None)]
    for atom in to_add:
        assert isinstance(atom, Atom), "to_add has invalid: %s" % str(atom)
        add = _get_next_subsmiles(atom, start_atom, molecule)
        output += "("+BOND_SYMBOLS[start_atom.non_h_neighbors[atom]]+add+")"

    return output



def _initialize_nonHNeighbors(molecule):
    "Create the dictionary nonHNeighbors for each atom."
    for atom in molecule.atoms:
        # atom.non_h_neighbors=dict((a,b) for (a,b)
        #       in atom.neighbors.items() if a.element.lower()!="h")
        ## TODO? Ignoring hydrogens is lame.
        atom.non_h_neighbors = copy.copy(atom.neighbors)

def _flag_rings(molecule):
    "Traverse the molecule once, to hunt down and flag rings."
    ringsfound = 0
    atom = molecule.atoms[0] ## current atom
    home = molecule.atoms[0] ## home atom

    #Each iteration: (...while we aren't back to the home atom, or if we are,
                    #while the home atom still has neighbors to read)
    while atom != home or home.n_read < len(home.non_h_neighbors):
    
        #flag current atom as "read" (flag = 1)
        atom.flag = 1
        
        #if there are neighbors left to read from this atom:
        if atom.n_read < len(atom.non_h_neighbors):
            
            neighbors = list(atom.non_h_neighbors)

            #if the next atom is the parent atom:
            if neighbors[atom.n_read] == atom.parent_atom:
                #don't do anything but incrementing nRead
                atom.n_read += 1
                
            #else,
            else:
                if neighbors[atom.n_read].n_read == 0:
                #if the next atom has not been traversed already:
                    atom.n_read += 1
                    neighbors[atom.n_read - 1].parent_atom = atom
                    atom = neighbors[atom.n_read - 1]
                    #increment current atom's nRead counter
                    #make the next atom the current atom:
                    #make the old atom the next atom's parent

                else:
                #if the next atom has been traversed already:
                    #it's a ring!
                    ringsfound += 1
                    atom.rflag.append((
                        ringsfound,
                        neighbors[atom.n_read]
                    ))
                    neighbors[atom.n_read].rflag.append((
                        ringsfound, atom
                    ))
                    atom.n_read += 1
                    #increment ringsfound
                    #set rflag on both atoms to ringsfound
                    #increment current atom's nRead counter
                    #make sure the atoms know who each other is
                
        #if not:
        else:
            #go backwards to parent atom:
            #set atom to its parent atom
            atom = atom.parent_atom

def _get_generated_smiles(molecule):
    """Precondition: Molecule has been flagged for rings already,
    using _flag_rings. Traverses a second time to generate SMILES."""
    start_atom = molecule.atoms[0]
    ## TODO: This is for debugging, and is temporary
    for atom in molecule.atoms:
        if atom.element == 'N':
            start_atom = atom
    return _subsmiles(molecule, start_atom, 0)

def _resetflags(molecule):
    "Reset all old flags."
    for atom in molecule.atoms:
        atom.flag = 0
        atom.rflag = []
        atom.n_read = 0
        atom.parent_atom = 0
        atom.non_h_neighbors = []


def _rflag_to_str(rflag):
    """
    Convert 3 to 3 and 45 to %45. For ring bonds.
    rflag :: int
    return :: str
    """
    assert isinstance(rflag, int), "Not an int rflag: %s" % str(rflag)
    if 0 <= rflag <= 9:
        return str(rflag)
    elif rflag <= 99:
        return "%" + str(rflag)
    else:
        ## The SMILES spec requires this. We could do it better by recycling
        ## ring-flags intelligently. But that would require a lot of work. TODO.
        raise StandardError("Too many rings in molecule. 100 is too many.")


def _get_subsmiles_ct(output, molecule, start_atom, parent_atom):
    """
    Check if the atom is a cis-trans center. Output correctly if so.
    Remember to worry about cis-trans centers that might be in a ring system.
    Remember to worry about whether or not an atom has a parent atom.
    Adds ring labels.
    """
    atomsToLink = [start_atom.CTotherC, start_atom.CTa, start_atom.CTb]
    begin = ["", "/", "\\"]
    if start_atom.CTotherC.flag == 2:
        begin = ["", "\\", "/"]
    if start_atom.CTa == parent_atom:
        output = begin[2] + output
    if start_atom.CTb == parent_atom:
        output = begin[1] + output
    for ind in range(3):
        atom = atomsToLink[ind]
        if (atom != None) and (atom != parent_atom):
            if atom in [rf[1] for rf in start_atom.rflag]:
                output += ''.join([
                    "(",
                    begin[ind],
                    BOND_SYMBOLS[start_atom.non_h_neighbors[atom]],
                    _rflag_to_str(start_atom.rflag[
                        [rf[1] for rf in start_atom.rflag].index(atom)][0]),
                    ")",
                ])
            elif atom.flag == 1:
                output += ''.join([
                    "(",
                    begin[ind],
                    BOND_SYMBOLS[start_atom.non_h_neighbors[atom]],
                    _subsmiles(molecule, atom, start_atom),
                    ")",
                ])
    return output


def _get_next_subsmiles(atom, start_atom, molecule):
    """
    Recursion is your friend.
    Be sure to specify the base case (when zero non-parent non-ring atoms are 
        available to bond to).
    In the base case, this loop won't even be entered.
    """
    if (start_atom.rflag != []) and (atom in [r[1] for r in start_atom.rflag]):
        return _rflag_to_str(start_atom.rflag[
            [r[1] for r in start_atom.rflag].index(atom)][0])
    else:
        return _subsmiles(molecule, atom, start_atom)


def _assertMolecule(molecule):
    "Assert that molecule is a Molecule."
    try:
        assert isinstance(molecule, Molecule)
    except AssertionError:
        raise StandardError("Not a molecule: %s is a %s, not %s" % \
            (repr(molecule), type(molecule), Molecule))
