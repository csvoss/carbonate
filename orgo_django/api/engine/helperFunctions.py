from molecularStructure import *
from toSmiles import *
import copy
import itertools
import cPickle

randThing = 0

debug = False


#Returns a list of molecules.
def antiAdd(molecule, target1, target2, add1, add2,
            addtarget1 = None, addtarget2 = None):
    #Just a wrapper function, to make antiAdd less confusing.
    return synAdd(molecule, target1, target2, add1, add2,
                  addtarget1, addtarget2, True)

def duplicateInputs(molecule, target1, target2, add1, add2, addtarget1,
                    addtarget2):
    #Helper function for adds.  Returns a deep-copied version of all the inputs.
    target1Pos = molecule.atoms.index(target1)
    target2Pos = molecule.atoms.index(target2)
    if addtarget1 != None:
        addtarget1Pos = add1.atoms.index(addtarget1)
    if addtarget2 != None:
        addtarget2Pos = add2.atoms.index(addtarget2)

    Xmolecule = copy.deepcopy(molecule)
    #Remake pointers to targets
    Xtarget1 = Xmolecule.atoms[target1Pos]
    Xtarget2 = Xmolecule.atoms[target2Pos]
    Xadd1 = copy.deepcopy(add1)
    Xadd2 = copy.deepcopy(add2)
    if addtarget1 != None:
        Xaddtarget1 = Xadd1.atoms[addtarget1Pos]
    else:
        Xaddtarget1 = None
    if addtarget2 != None:
        Xaddtarget2 = Xadd2.atoms[addtarget2Pos]
    else:
        Xaddtarget2 = None
    return (Xmolecule, Xtarget1, Xtarget2, Xadd1, Xadd2, Xaddtarget1, Xaddtarget2)

#Returns a list of molecules.
def synAdd(molecule, target1, target2, add1, add2,
           addtarget1 = None, addtarget2 = None, antiAdd = False):
    #Destroys the double bond and CTstereochemistry between target1 and target2.
    #Adds add1 and add2 to target1 and target2.  If add1 and/or add2 are molecules,
    #addtargets are needed to specify where the bond should originate from add.

    #Also does anti-addition, if antiAdd is set to true.
    (molecule, target1, target2, add1, add2, addtarget1, addtarget2) =\
               duplicateInputs(molecule, target1, target2, add1, add2, addtarget1,
                               addtarget2)
    (Xmolecule, Xtarget1, Xtarget2, Xadd1, Xadd2, Xaddtarget1, Xaddtarget2) =\
               duplicateInputs(molecule, target1, target2, add1, add2, addtarget1,
                               addtarget2)

    #Set bond orders to single.
    target1.neighbors[target2] = 1
    target2.neighbors[target1] = 1
    Xtarget1.neighbors[Xtarget2] = 1
    Xtarget2.neighbors[Xtarget1] = 1

    if antiAdd:
        bigListOfStuff =\
        ((molecule, add1, target1, addtarget1, target2, target1.CTa, target1.CTb),
        (molecule, add2, target2, addtarget2, target1, target2.CTa, target2.CTb),
        (Xmolecule, Xadd1, Xtarget1, Xaddtarget1, Xtarget2, Xtarget1.CTb, Xtarget1.CTa),
        (Xmolecule, Xadd2, Xtarget2, Xaddtarget2, Xtarget1, Xtarget2.CTb, Xtarget2.CTa))
    else:
        bigListOfStuff =\
        ((molecule, add1, target1, addtarget1, target2, target1.CTb, target1.CTa),
        (molecule, add2, target2, addtarget2, target1, target2.CTa, target2.CTb),
        (Xmolecule, Xadd1, Xtarget1, Xaddtarget1, Xtarget2, Xtarget1.CTa, Xtarget1.CTb),
        (Xmolecule, Xadd2, Xtarget2, Xaddtarget2, Xtarget1, Xtarget2.CTb, Xtarget2.CTa))

    for thismolecule, thisAdd, thisTarget, thisAddTarget, otherTarget, ct1, ct2\
            in bigListOfStuff:
        if isinstance(thisAdd, Atom):
            thismolecule.addAtom(thisAdd, thisTarget, 1)
            if ct1 != None or ct2 != None:
                thisTarget.newChiralCenter(otherTarget,
                        (thisAdd, ct1, ct2))
        elif isinstance(thisAdd, Molecule):
            #Untested.
            thismolecule.addMolecule(thisAdd, thisAddTarget, thisTarget, 1)
            if ct1 != None or ct2 != None:
                thisTarget.newChiralCenter(otherTarget,
                        (thisAddTarget, ct1, ct2))
        else:
            #Hydrogens.
            if ct1 != None and ct2 != None:
                thisTarget.newChiralCenter(otherTarget,
                        (None, ct1, ct2))
        thisTarget.eliminateCT()
    if moleculeCompare(molecule, Xmolecule):
        return [molecule]
    else:
        return [molecule, Xmolecule]

def allAdd(molecule, target1, target2, add1, add2, addtarget1=None, addtarget2=None):
    #Adds add1 and add2 to target1 and target2, in both syn and anti fashions.
    #Does not introduce stereochemistry (this kind of addition never results in
    #stereochemistry).

    #Protect the inputs from modification:
    (molecule, target1, target2, add1, add2, addtarget1, addtarget2)=\
               duplicateInputs(molecule, target1, target2, add1, add2, addtarget1, addtarget2)
    #Reduces bond order by 1
    bo = target1.neighbors[target2]
    molecule.changeBond(target1, target2, bo-1)
    target1.eliminateCT()
    target2.eliminateCT()
    target1.eliminateChiral()
    target2.eliminateChiral()
    #Add new stuff
    for thisTarget, thisAdd, thisAddtarget in (
        (target1, add1, addtarget1), (target2, add2, addtarget2)):
        if isinstance(thisAdd, Molecule):
            molecule.addMolecule(thisAdd, thisAddtarget, thisTarget, 1)
        elif isinstance(thisAdd, Atom):
            molecule.addAtom(thisAdd, thisTarget, 1)
        else:
            #Hydrogens.  Do nothing.
            pass
    return [molecule]

def tripleAdd(molecule, target1, target2, add1, add2, cisOrTrans,
              addtarget1=None, addtarget2=None):
    #A triple bond must be between targets 1 and 2.  Adds add1 and add2
    #while reducing order by 1 (to double bond).  cisOrTrans must be specified,
    #as a string ("cis" or "trans")
    #As always, addtargets are needed iff add1 and add2 are molecules.
    
    #Protect the inputs from modification:
    (molecule, target1, target2, add1, add2, addtarget1, addtarget2)=\
        duplicateInputs(molecule, target1, target2, add1, add2, addtarget1, addtarget2)

    #Change bond orders
    if target1.neighbors[target2] != 3:
        print "Error in tripleAdd: no triple bond specified."
        raise StandardError

    stuff = ((molecule, target1, target2, add1, addtarget1),
     (molecule, target2, target1, add2, addtarget2))
    for thisMol, thisTarget, otherTarget, thisAdd, thisAddtarget in stuff:
        thisMol.changeBond(thisTarget, otherTarget, 2)
        if isinstance(thisAdd, Atom):
            thisMol.addAtom(thisAdd, thisTarget, 1)
            CTthing = thisAdd
        elif isinstance(thisAdd, Molecule):
            thisMol.addMolecule(thisAdd, thisAddtarget, thisTarget, 1)
            CTthing = thisAddtarget
        else:
            #Hydrogen
            CTthing = None
        otherAttached = None #By default, otherAttached is a hydrogen.
        for neighbor in thisTarget.neighbors:
            if neighbor != otherTarget and neighbor != thisAdd and neighbor != thisAddtarget:
                otherAttached = neighbor
        if cisOrTrans.lower() == 'trans':
            thisTarget.newCTCenter(otherTarget, otherAttached, CTthing)
        else:
            if hasattr(otherTarget, "CTotherC"):
                thisTarget.newCTCenter(otherTarget, otherAttached, CTthing)
            else:
                thisTarget.newCTCenter(otherTarget, CTthing, otherAttached)
    return [molecule]
 
def allTripleAdd(molecule, target1, target2, add1, add2, addtarget1 = None, addtarget2 = None):
    #Adds two copies of add1 and two copies of add2 to target1 and target2, respectively.
    #Breaks a triple bond.  Introduces no new stereochemistry.

    #Protect the inputs from modification:
    (molecule, target1, target2, add1, add2, addtarget1, addtarget2)=\
               duplicateInputs(molecule, target1, target2, add1, add2, addtarget1, addtarget2)
    #We need an extra copy of add1 and add2, along with corresponding addtargets.
    (notused, notused2, notused3, add1b, add2b, addtarget1b, addtarget2b)=\
               duplicateInputs(molecule, target1, target2, add1, add2, addtarget1, addtarget2)
    #Change to single bond
    molecule.changeBond(target1, target2, 1)
    #Add new stuff
    for thisTarget, thisAdd, thisAddtarget in (
        (target1, add1, addtarget1), (target2, add2, addtarget2),
        (target1, add1b, addtarget1b), (target2, add2b, addtarget2b)):
        if isinstance(thisAdd, Molecule):
            molecule.addMolecule(thisAdd, thisAddtarget, thisTarget, 1)
        elif isinstance(thisAdd, Atom):
            molecule.addAtom(thisAdd, thisTarget, 1)
        else:
            #Hydrogens.  Do nothing.
            pass
    return [molecule]
    
def moleculeCompare(a, b, checkChiral = True):
    global finishedDict
    finishedDict = None
    def insideCompare(a, b, compareDict = None, expanded = []):
        global finishedDict
        #Determines whether two molecules are isomorphic.  In the worst case
        #(two molecules with the same atoms), this procedure does not run in
        #polynomial time, so be careful.
        for ele in ['C','N','O']:
            if a.countElement(ele) != b.countElement(ele):
                return False
        #sa = smiles(a)
        #sb = smiles(b)

        #compareDict maps atoms in a to their hypothesized counterparts in b.

        if len(expanded) == len(a.atoms):
            #We've reached every atom.  Call it equal.
            finishedDict = compareDict
            return True

        if compareDict == None:
            for atom in b.atoms:
                if atom.element == a.atoms[0].element:
                    oldCompareDict = {a.atoms[0]:atom, None:None}
                    newCompareDicts = neighborCompare(a.atoms[0], atom,
                                    oldCompareDict)
                    if newCompareDicts == None:
                        continue
                    for newCompareDict in newCompareDicts:
                        if insideCompare(a, b, dict(newCompareDict.items()+
                                        oldCompareDict.items()), [a.atoms[0]]):
                            return True
            return False

        for aAtom in compareDict:
            if (aAtom in expanded) or aAtom == None:
                #Already expanded this atom.  Don't do it again.
                continue
            newDictSectors = neighborCompare(aAtom, compareDict[aAtom], compareDict)
            if newDictSectors == None:
                return False
            for newDictSector in newDictSectors:
                if insideCompare(a, b, dict(compareDict.items() + newDictSector.items()),
                                   expanded + [aAtom]):
                    return True
            return False
                
    def neighborCompare(a,b, compareDict):
        #Helper function.  Given 2 atoms, returns all pairings of neighbors of a
        #with neighbors of b such that each pair has the same element.
        aN = []
        bN = []
        for aNeighbor in a.neighbors:
            aN.append(aNeighbor.element)
        for bNeighbor in b.neighbors:
            bN.append(bNeighbor.element)
        #If the elements don't match, obviously there are no pairings.
        if sorted(aN) != sorted(bN):
            return None
        if hasattr(a, "chiralA") != hasattr(b, "chiralA"):
            #One atom has chirality, where the other doesn't.  Obviously no pairings.
            return None
        if hasattr(a, "chiralA"):
            chiralFlag = True
        else:
            chiralFlag = False
        if hasattr(a, "CTotherC") != hasattr(b, "CTotherC"):
            #One atom has a cis-trans center, where the other doesn't.  Obviously no pairings.
            return None
        if hasattr(a, "CTotherC"):
            CTFlag = True
        else:
            CTFlag = False
        #Generate all n! pairings, and prune as we go.
        out = []
        for aNeighborSet in itertools.permutations(a.neighbors):
            temp = {None: None}
            OKFlag = True
            for i in xrange(len(aNeighborSet)):
                if aNeighborSet[i].element != b.neighbors.keys()[i].element:
                    #Oops, the elements don't actually match.  Skip.
                    OKFlag = False
                    break
                if a.neighbors[aNeighborSet[i]] != b.neighbors[b.neighbors.keys()[i]]:
                    #Oops, the bond orders are different.  Skip.
                    OKFlag = False
                    break
                if aNeighborSet[i] in compareDict and \
                   b.neighbors.keys()[i] != compareDict[aNeighborSet[i]]:
                    #This pairing goes against what's already in compareDict.  Skip.
                    OKFlag = False
                    break
                temp[aNeighborSet[i]] = b.neighbors.keys()[i]
            if chiralFlag and OKFlag and checkChiral:
    #            if not(set(a.neighbors.keys()) <= set((a.chiralA, a.chiralB, a.chiralC, a.chiralD))):
    #                if debug:
    #                    print "Error thing"
    #                    for neighbor in a.neighbors.keys():
    #                        if neighbor != None:
    #                            print neighbor.element
    #                    print "------"
    #                    for ch in (a.chiralA, a.chiralB, a.chiralC, a.chiralD):
    #                        if ch != None:
    #                            print ch.element
    #                    print "------"
    #                    raise StandardError
                #The following bit of code is still quite messy.  It tests whether the
                #hypothesized pairing follows the correct chirality.
                aCW = []
                bCW = []
                for neighbor in a.chiralCWlist(aNeighborSet[randThing]):
                    if neighbor == None:
                        #Hydrogen.
                        aCW.append("H")
                    else:
                        aCW.append(neighbor.element)
                for neighbor in b.chiralCWlist(b.neighbors.keys()[randThing]):
                    if neighbor == None:
                        #Hydrogen.
                        bCW.append("H")
                    else:
                        bCW.append(neighbor.element)
                OKFlag = False
                for i in xrange(3):
                    #Find the correct alignment of neighbors of a to neighbors of b
                    if OKFlag:
                        break
                    if aCW == shift(bCW, i):
                        OKFlag = True
                        #The elements are correct, but is the actual mapping consistant?
                        for j in xrange(3):
                            if a.chiralCWlist(aNeighborSet[randThing])[j] == None:
                                continue
                            if temp[a.chiralCWlist(aNeighborSet[randThing])[j]] ==\
                                b.chiralCWlist(b.neighbors.keys()[randThing])[(j+i)%3]:
                                pass
                            else:
                                OKFlag = False
            if CTFlag and OKFlag:
                #Makes sure that the hypothesized pairing follows the correct
                #cis-trans relationship
                if a.CTotherC in compareDict and a.CTotherC.CTa in compareDict and a.CTotherC.CTb in compareDict:                    
                    if ((b.CTa == temp[a.CTa]) !=
                       (b.CTotherC.CTa == compareDict[a.CTotherC.CTa])) or\
                       ((b.CTb == temp[a.CTb]) !=
                        (b.CTotherC.CTb == compareDict[a.CTotherC.CTb])):
                        OKFlag = False
                    #If any C has 2 hydrogens, the matching is automatically correct, as there is no stereochem.
                    if (a.CTa == None and a.CTb == None) or (a.CTotherC.CTa == None and a.CTotherC.CTb == None):
                        OKFlag = True

            del(temp[None])
            if (temp not in out) and OKFlag:
                out.append(temp)
        return out
        
    if checkChiral:
        return insideCompare(a, b)
    else:
        ans = insideCompare(a, b)
        if (finishedDict == None) and ans:
            raise StandardError
        return ans, finishedDict

def shift(l, n):
    return l[n:] + l[:n]

def isInRing(start, current, last=None, alreadyVisited=[]):
    #Tests whether ADJACENT Atoms start and current are in the same ring.
    #Only works on adjacent atoms!
    #If so, returns a list of elements in the ring in connected order, starting (current, ..., start)
    #If not, returns None
    #Do a depth-first search.
    if last == None:
        #Initialization only.
        last = start
    if start == current:
        return [current]
    if current in alreadyVisited:
        return None
    for nextNeighbor in current.neighbors:
        if nextNeighbor == last:
            #Don't go where we've already been.
            continue
        x = isInRing(start, nextNeighbor, current, alreadyVisited + [current])
        if x!= None:
            return [current] + x
    return None

def markovnikov(a, b):
    #a and b are two carbon atoms.  Function tuple of all possible markovnikov
    #orderings of carbons
    aTotal = 0
    bTotal = 0
    for atom in a.neighbors:
        if atom.element == "C":
            aTotal += 1
    for atom in b.neighbors:
        if atom.element == "C":
            bTotal += 1
    if aTotal == bTotal:
        return ((a,b),(b,a))
    elif aTotal > bTotal:
        return ((a, b),)
    else:
        return ((b, a),)

#Converts an alkyne to a carbonyl at target1 carbon.
def carbonylAdd(molecule, target1, target2):
    add1 = Atom("O")
    add2 = None
    addtarget1 = None
    addtarget2 = None

    #Protect the inputs from modification:
    (molecule, target1, target2, add1, add2, addtarget1, addtarget2)=\
               duplicateInputs(molecule, target1, target2, add1, add2, addtarget1, addtarget2)
    #We need an extra copy of add1 and add2, along with corresponding addtargets.
    (notused, notused2, notused3, add1b, add2b, addtarget1b, addtarget2b)=\
               duplicateInputs(molecule, target1, target2, add1, add2, addtarget1, addtarget2)
    #Change to single bond
    molecule.changeBond(target1, target2, 1)
    #Add new stuff
    molecule.addAtom(add1, target1, 2)
    
    return [molecule]


#Returns a tuple of atoms.
#Returns None if none found.
def findAlkene(molecule):
    for atom in molecule.atoms:
        if not (atom.element == 'C'):
            continue
        for neighbor in atom.neighbors:
            if neighbor.element == 'C' and atom.neighbors[neighbor] == 2:
                return (atom, neighbor)
    return None


#Returns a list of tuples of atoms.
#Returns [] if none found.
def findAlkenes(molecule):
    output = []
    for atom in molecule.atoms:
        if not (atom.element == 'C'):
            continue
        for neighbor in atom.neighbors:
            if neighbor.element == 'C' and atom.neighbors[neighbor] == 2:
                output += [(atom, neighbor)]
    #Remove duplicate tuples
    output2 = []
    for item in output:
        if not ((item[1], item[0]) in output2):
            output2 += [item]
    return output2



#Returns a tuple of atoms.
#Returns None if none found.
def findAlkyne(molecule):
    if molecule == None:
        return None
    for atom in molecule.atoms:
        if not (atom.element == 'C'):
            continue
        for neighbor in atom.neighbors:
            if neighbor.element == 'C' and atom.neighbors[neighbor] == 3:
                if atom == neighbor:
                    print "OH NO WHAT JUST HAPPENED"
                    raise StandardError
                return (atom, neighbor)
    return None


#Returns an atom. That atom is part of an alkyne, and is attached to an H.
#Returns None if none found.
def findHydrogenAlkyne(molecule):
    for atom in molecule.atoms:
        if not (atom.element == 'C'):
            continue
        if len(list(atom.neighbors)) != 1:
            continue
        if atom.charge != 0:
            continue
        for neighbor in atom.neighbors:
            if neighbor.element == 'C' and atom.neighbors[neighbor] == 3 and neighbor.charge == 0:
                return atom
    return None

#Returns a list of tuples of atoms.
#Returns [] if none found.
def findAlkynes(molecule):
    output = []
    for atom in molecule.atoms:
        if not (atom.element == 'C'):
            continue
        for neighbor in atom.neighbors:
            if neighbor.element == 'C' and atom.neighbors[neighbor] == 3:
                output += [(atom, neighbor)]
    #Remove duplicate tuples
    output2 = []
    for item in output:
        if not ((item[1], item[0]) in output2):
            output2 += [item]
    return output2
    
#Returns a list of atoms.
#Returns [] if none found.
def findHydroxyls(molecule):
    output = []
    if moleculeCompare(Molecule(Atom("O")), molecule):
        return [[x for x in molecule.atoms if x.element == "O"][0]]
    for atom in molecule.atoms:
        if not (atom.element == 'O'):
            continue
        if len(list(atom.neighbors)) == 0: #water
            output += [atom]
        elif len(list(atom.neighbors)) == 1: #it's a hydroxyl?
            if atom.neighbors.values()[0] == 1: #single bonds only!  No ketones!
                output += [atom]
            else:
                continue
        elif len(list(atom.neighbors)) == 2: #it's an ether, don't return it
            continue
        else:
            print "Error -- Invalid oxygen atom with 3+ neighbors."
            raise StandardError
    return output

#Returns a list of tuples of atoms.
#Returns [] if none found.
def findAlkenesAndAlkynes(molecule):
    return findAlkenes(molecule) + findAlkynes(molecule)

#Returns a tuple of atoms.
#Returns None if none found.
def findAlkeneAndAlkyne(molecule):
    #Tiny helper function.
    x = findAlkene(molecule)
    if x == None:
        return findAlkyne(molecule)
    else:
        return x

def findAlkeneOrAlkyne(molecule):
    return findAlkeneAndAlkyne

#def findAlkenes(molecule):
#    return findAlkenesOrAlkynes(molecule, 2)

#def findAlkynes(molecule):
#    return findAlkenesOrAlkynes(molecule, 3)

#Finds all candidate alkenes within a molecule.
#(define "alkenes" as "alkenes that are not in an aromatic ring")
#Returns a list of tuples of atoms. The lowest tuple is a pair of two atoms, which share a double bond.
#Make sure not to include duplicates.
def findAlkenesOrAlkynes(molecule, bo):
    #bo = bond order (2 or 3)
    #To track which bonds we've counted, we use the atom.flag property.
    #atom.flag starts at 0, and must be reset to 0 at the end.
    doubleBonds = []
    for atom in molecule.atoms:
        if not(atom.element == 'C'):
            continue
        for neighbor in atom.neighbors:
            if neighbor.element == 'C' and atom.neighbors[neighbor] == bo:
                if atom.flag != 0 and neighbor.flag != 0 and atom in neighbor.flag\
                   and neighbor in atom.flag:
                    #We've already counted this bond.  Move on.
                    continue
                if atom.flag == 0:
                    atom.flag = [neighbor]
                else:
                    atom.flag.append(neighbor)
                if neighbor.flag == 0:
                    neighbor.flag = [atom]
                else:
                    neighbor.flag.append(atom)
                doubleBonds.append((atom, neighbor))
    #Reset the flags - other functions like to use them, as well.
    for atom in molecule.atoms:
        atom.flag = 0
    return doubleBonds








#Takes in a molecule or list of molecules
#Outputs a list of molecule(s), representing distinct contiguous parts of the old molecule(s).
#Uses include ozonolysis reactions.
def splice(molecules):
    if not isinstance(molecules, list):
        molecules = [molecules]
    output = []
    for molecule in molecules:
        #Find the next unflagged atom, if it exists
        while 0 in [atom.flag for atom in molecule.atoms]:
            originalAtom = [atom for atom in molecule.atoms if atom.flag == 0][0]
            currentAtom = originalAtom
            newMolecule = Molecule(originalAtom)
            newMolecule.atoms = []
            #Traverse the map of bonds, storing parent atoms and checking for flags
            #Do this until you are back at the original atom with nothing left to check
            while (0 in [atom.flag for atom in list(originalAtom.neighbors)]) or (currentAtom != originalAtom):
                #If current atom unflagged:
                if currentAtom.flag == 0:
                    #Add current atom to list
                    newMolecule.atoms += [currentAtom]
                    #Flag current atom
                    currentAtom.flag = 1
                #If current atom lacks unflagged neighbors:
                if not 0 in [atom.flag for atom in list(currentAtom.neighbors)]:
                    #Return to parent atom.
                    currentAtom = currentAtom.parentAtom
                #Else:
                else:
                    #Next atom is the first unflagged neighbor.
                    nextAtom = [atom for atom in list(currentAtom.neighbors) if atom.flag == 0][0]
                    #Store current atom as next atom's parent.
                    nextAtom.parentAtom = currentAtom
                    #Store next atom as current atom.
                    currentAtom = nextAtom
            output += [newMolecule]


        #Reset flags
        for atom in molecule.atoms:
            atom.flag = 0
            atom.parentAtom = 0
    return output
    
    
def verify(molecule):
    #Checks to make sure that the linkage and stereochem of each atom match.
    for atom in molecule.atoms:
        if hasattr(atom, 'chiralA'):
            for neighbor in (atom.chiralA, atom.chiralB, atom.chiralC, atom.chiralD):
                if neighbor == None:
                    continue
                if neighbor not in atom.neighbors.keys():
                    print "-----Error - chirality broken A!-----"
            for neighbor in atom.neighbors.keys():
                if neighbor not in (atom.chiralA, atom.chiralB, atom.chiralC, atom.chiralD):
                    print "-----Error - chirality broken B!-----"
        if hasattr(atom, 'CTa'):
            for neighbor in (atom.CTotherC, atom.CTa, atom.CTb):
                if neighbor == None:
                    continue
                if neighbor not in atom.neighbors.keys():
                    print "-----Error - CT broken A!-----"
            for neighbor in atom.neighbors.keys():
                if neighbor not in (atom.CTotherC, atom.CTa, atom.CTb):
                    print "-----Error - CT broken B!-----"



