"""
toSmiles.py

This code contains the function smilesify(molecule), which converts a Molecule object to a canonical SMILES string:
http://en.wikipedia.org/wiki/SMILES
Once a Molecule has been converted to a SMILES string, we can render it as an SVG using renderSVG.py.

Also contains the function to_canonical(smiles), which converts a SMILES string to a canonical SMILES string.

PROCEED CAUTIOUSLY -- THIS CODE IS FINICKY

Does not yet support:
- charge :( :( 
- radicals :( :( 
- isotopes
"""
from molecularStructure import *

import openbabel
import pybel

debugSmiles = False


def to_canonical(smiles):
    """
    Uses OpenBabel.
    Converts a SMILES string to a canonical SMILES string.
    smiles :: str.
    return :: str.
    """
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("smi", "can")
    outMol = openbabel.OBMol()
    obConversion.ReadString(outMol, str(smiles))
    ans = obConversion.WriteString(outMol)
    return ans

def assertMolecule(molecule):
    try:
        assert isinstance(molecule, Molecule)
    except AssertionError:
        raise StandardError("Not a molecule: %s is a %s, not %s" % (repr(molecule), type(molecule), Molecule) )

def smilesify(molecule):
    """
    molecule :: Molecule or [Molecule]
    (in the latter case, smilesify(mol) is applied for mol in the list)

    Converts Molecule objects to SMILES strings. Traverses the molecule to detect and mark
    rings or cycles. Uses a convoluted system of flags to do this. Then, passes it on to 
    subsmiles, which operates on a molecule with rings already flagged, and performs tree
    traversal.

    return :: str or [str]
    """

    if isinstance(molecule, list):
        try:
            assertMolecule(molecule[0])
        except IndexError:
            return []
        return [smilesify(molec) for molec in molecule]



    if len(molecule.atoms)==0: return ""
    
    ringsfound = 0
    cur_atom = molecule.atoms[0]
    home_atom = molecule.atoms[0]

    #Create the dictionary nonHNeighbors for each atom.
    for atom in molecule.atoms:
        atom.nonHNeighbors = dict((a, b) for (a, b)in atom.neighbors.items() if a.element.lower() != "h")
    

    #Traverse the molecule once, to hunt down and flag rings.
    #Each iteration: (...while we aren't back to the home atom, or if we are,
                    #while the home atom still has neighbors to read)
    while ((cur_atom != home_atom) or (home_atom.nRead < len(home_atom.nonHNeighbors))):
    
        #flag current atom as "read" (flag = 1)
        cur_atom.flag = 1
        if debugSmiles:
            print "Flagged "+str(cur_atom)+" as read."
        
        #if there are neighbors left to read from this atom:
        if (cur_atom.nRead < len(cur_atom.nonHNeighbors)):
            
            #if the next atom is the parent atom:
            if list(cur_atom.nonHNeighbors)[cur_atom.nRead] == cur_atom.parentAtom:
                #don't do anything but incrementing nRead
                cur_atom.nRead += 1
                if debugSmiles:
                    print "Nope, not progressing to parent."
                
            #else,
            else:
                if list(cur_atom.nonHNeighbors)[cur_atom.nRead].nRead == 0:
                #if the next atom has not been traversed already:
                    cur_atom.nRead += 1
                    list(cur_atom.nonHNeighbors)[cur_atom.nRead - 1].parentAtom = cur_atom
                    cur_atom = list(cur_atom.nonHNeighbors)[cur_atom.nRead - 1]
                    #increment current atom's nRead counter
                    #make the next atom the current atom:
                    #make the old atom the next atom's parent
                    if debugSmiles:
                        print "Progressing..."

                else:
                #if the next atom has been traversed already:
                    #it's a ring!
                    ringsfound += 1
                    cur_atom.rflag += [(ringsfound, list(cur_atom.nonHNeighbors)[cur_atom.nRead])]
                    list(cur_atom.nonHNeighbors)[cur_atom.nRead].rflag += [(ringsfound, cur_atom)]
                    cur_atom.nRead += 1
                    #increment ringsfound
                    #set rflag on both atoms to ringsfound
                    #increment current atom's nRead counter
                    #make sure the atoms know who each other is
                    if debugSmiles:
                        print "Ring "+str(ringsfound)+" found! "+str(cur_atom.rflag)+", "+str(list(cur_atom.nonHNeighbors)[cur_atom.nRead - 1])
                
        #if not:
            #go backwards to parent atom:
            #set cur_atom to its parent atom
        else:
            if debugSmiles:
                print "Regressing to "+str(cur_atom.parentAtom)+" from "+str(cur_atom)
            cur_atom = cur_atom.parentAtom
            
            
    #Traverse twice to generate the SMILES.
            
        
    if debugSmiles:
        for atom in molecule.atoms:
            print ""
            print atom.element
            print atom
            print atom.rflag
            print atom.nonHNeighbors
        
    outp = subsmiles(molecule, molecule.atoms[0], 0)

    #Reset all old flags.    
    for atom in molecule.atoms:
        atom.flag = 0
        atom.rflag = []
        atom.nRead = 0
        atom.parentAtom = 0
        atom.nonHNeighbors = []

    return to_canonical(outp)

bondSymbols = ['0', '-', '=', '#', '$']

def subsmiles(molecule, startAtom, parentAtom):

    """
    Precondition: molecule has been flagged for ring positioning (some rflag values on atoms might != 0). This is done by smilesify().
    Creates and returns a SMILES string for unflagged (!atom.flag==2) atoms within a molecule, starting with the given atom.

    molecule :: Molecule.
    startAtom :: Atom.
    parentAtom :: Atom.

    Traverses the molecule from the given starting atom, returning the SMILES representation.
    Can be called recursively -- hence startAtom and parentAtom. This is a tree traversal, after all!

    return :: a SMILES substring
    """

    if debugSmiles:
        if parentAtom != 0:
            print "\t\t\t [ Entering subsmiles with "+str((startAtom.element, startAtom, startAtom.rflag))+" from "+str((parentAtom.element, parentAtom, parentAtom.rflag))
        else:
            print "\t\t\t [ Entering subsmiles with "+str((startAtom.element, startAtom, startAtom.rflag))+"..."
    
    #Flag the current atom.
    startAtom.flag = 2
    
    outp = str(startAtom.element)

    #Check if the atom is a cis-trans center. Output correctly if so.
    #Remember to worry about cis-trans centers that might be part of a ring system.
    #Remember to worry about whether or not an atom has a parent atom.
    #Adds ring labels.
    if hasattr(startAtom, 'CTotherC'):
        if debugSmiles:
            print "\t\t\tSubsmiles case - CT"
        #print (True, startAtom.element)
        atomsToLink = [startAtom.CTotherC, startAtom.CTa, startAtom.CTb]
        begin = ["", "/", "\\"]
        if startAtom.CTotherC.flag == 2:
            begin = ["", "\\", "/"]
        if startAtom.CTa == parentAtom:
            outp = begin[2] + outp
        if startAtom.CTb == parentAtom:
            outp = begin[1] + outp
        for ind in range(3):
            atom = atomsToLink[ind]
            if (atom != None) and (atom != parentAtom):
                if atom in [rf[1] for rf in startAtom.rflag]:
                    outp += "(" + begin[ind] + bondSymbols[startAtom.nonHNeighbors[atom]] + str(startAtom.rflag[[rf[1] for rf in startAtom.rflag].index(atom)][0]) + ")"
                elif atom.flag == 1:
                    outp += "(" + begin[ind] + bondSymbols[startAtom.nonHNeighbors[atom]] + subsmiles(molecule, atom, startAtom) + ")"
        if debugSmiles:
            print "\t\t\t >>> " + outp + " ] "
        return outp
    
   

    ###Put a ring marker on the atom, if its ring partner is not flagged yet.
    ##if (startAtom.rflag != 0) and (startAtom.rAtom.flag != 2):
    ##    outp += str(startAtom.rflag)

    #Check if the atom is a chiral center. If so:
    if hasattr(startAtom, 'chiralA'):
        if debugSmiles:
            print "\t\t\tSubsmiles case -- chiral"
        hasP = (parentAtom != 0)
        hasH = (True in [a.element.lower()=="h" for a in list(startAtom.neighbors)]) or (None in [startAtom.chiralA, startAtom.chiralB, startAtom.chiralC, startAtom.chiralD])
        #If the atom has a hydrogen:
        #Add [ and @@H] to the current output. (e.g. [C@@H]
        #If the atom does not have a hydrogen:
        #Add [ and @@] to the current output. (e.g. [C@@]
        #Then, add the remaining neighbors in the proper order.
        #Implement this by smart rearrangement of toAdd.
        #Be mindful of: whether or not there is a parent atom; whether or not there is a hydrogen.
        if hasH:
            outp = "[" + outp + "@@H]"
            if hasP:
                if debugSmiles:
                    print "\t\t\t...hasH hasP"
                #toAdd should have two elements
                l = startAtom.chiralCWlist(parentAtom) #list of three atoms
                x = l.index(None) #index of hydrogen atom in list
                toAdd = [l[(x+1) %3], l[(x+2) %3]] #correct permutation
                if None in toAdd:
                    print "\t\t\tError: atom "+startAtom.element+" is chiral, but has two hydrogens."
                    raise StandardError
            else:
                if debugSmiles:
                    print "\t\t\t...hasH !hasP"
                #toAdd should have three elements
                toAdd = startAtom.chiralCWlist(None) #list of three atoms
                if None in toAdd:
                    print "\t\t\tError: atom "+startAtom.element+" is chiral, but has two hydrogens."
                    raise StandardError
        else:
            outp = "[" + outp + "@@]"
            if hasP:
                if debugSmiles:
                    print "\t\t\t...!hasH hasP"
                #toAdd should have three elements
                toAdd = startAtom.chiralCWlist(parentAtom)
            else:
                if debugSmiles:
                    print "\t\t\t...!hasH !hasP"
                #toAdd should have four elements
                arbitraryRef = list(startAtom.neighbors)[0]
                l = startAtom.chiralCWlist(arbitraryRef)
                toAdd = [arbitraryRef] + l

    
    #Prepare to add new groups for all neighbor atoms which are not the parent atom and not the rAtom.
    else:
        if debugSmiles:
            print "\t\t\tSubsmiles case -- normal"
        toAdd = [atom for atom in list(startAtom.nonHNeighbors) if not (atom==parentAtom or atom==None)]

#    if toAdd == None:
#        toAdd = [atom for atom in list(startAtom.nonHNeighbors) if not (atom==parentAtom or atom==None)]

    
    #Recursion is your friend.
    #Be sure to specify the base case (when zero non-parent non-ring atoms are available to bond to)
    #In the base case, this loop won't even be entered.
    if debugSmiles:
        print "\t\t\t"+str([(atom.element, atom, atom.rflag) for atom in toAdd])
    for atom in toAdd:
        if (startAtom.rflag != []) and (atom in [rf[1] for rf in startAtom.rflag]) :
            add = str(startAtom.rflag[[rf[1] for rf in startAtom.rflag].index(atom)][0])
            if debugSmiles:
                print "Flagged as a ring bond -- "+startAtom.element+" to "+atom.element
        #elif atom.flag == 2:
            #You *really* don't want to enter subsmiles to this atom
            #pass
        else:
            if debugSmiles:
                print "Entering subsmiles, since "+str(startAtom.rflag != [])+" and "+str((atom in [rf[1] for rf in startAtom.rflag]))
            add = subsmiles(molecule, atom, startAtom)
                
        outp += "(" +bondSymbols[startAtom.nonHNeighbors[atom]] + add + ")"

    if debugSmiles:
        print "\t\t\t >>> " + outp + " ] "
    return outp