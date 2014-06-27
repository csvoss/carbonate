
from helperFunctions import *
import itertools

MAXLEN = 8

class ReactionTooCrazyError(Exception):
    pass

def removeDuplicates(moleculeList):
    if not isinstance(moleculeList, list):
        return [moleculeList]
    if len(moleculeList) == 0:
        return []
    return reduceChirality(removeDuplicatesAt(tautomerize(copy.deepcopy(moleculeList)), 0))


def removeDuplicatesAt(moleculeList, ind):
    if len(moleculeList) < ind+2:
        return moleculeList
    a = moleculeList[ind]
    for i in range(ind+1, len(moleculeList)):
        if a == moleculeList[i] or moleculeCompare(a, moleculeList[i]):
            del moleculeList[ind]
            return removeDuplicatesAt(moleculeList, ind)

    return removeDuplicatesAt(moleculeList, ind+1)

#Looks for two molecules that are the same, except at a single chiral center.
#Merges those two molecules by eliminating the chiral center.
#Assumes that molecules are all different.        
#MODIFIES THE INPUT molList
def reduceChirality(molList):
    if len(molList) == 1:
        return molList
    for a, b in itertools.combinations(molList, 2):
        same, compareDict = moleculeCompare(a, b, checkChiral=False)
        #In order for a center to be eliminated:
        #1) It has opposite chirality across a and b, while every other center has the same.
        nonMatchingCenter = None
        nonMatchCount = 0
        if same:
            #a and b are linked in the same way.
            for Acenter in compareDict.keys():
                Bcenter = compareDict[Acenter]
                if hasattr(Acenter, "chiralA") and hasattr(Bcenter, "chiralA"):
                    #Let's see if the chiralities are different.  We already know that the linkages
                    #around Acenter and Bcenter are the same.
                    Batoms1 = [compareDict[atom] for atom in Acenter.chiralCWlist(Acenter.chiralA)]
                    Batoms2 = Bcenter.chiralCWlist(compareDict[Acenter.chiralA])
                    #Are Batoms1 and Batoms2 equivalent up to a cycling?
                    if Batoms1 == Batoms2 or shift(Batoms1, 1) == Batoms2 or\
                        shift(Batoms1, 2) == Batoms2:
                        pass
                    else:
                        nonMatchingCenter = Acenter
                        nonMatchCount += 1
            if nonMatchCount == 1:
                #Remove the one non-matching center's chirality.
                molList.remove(b)
                nonMatchingCenter.eliminateChiral()
                return reduceChirality(molList)
        return molList


def react(molecules, findPlace, reactAtPlace):
    if not isinstance(molecules, list):
        return react([molecules], findPlace, reactAtPlace)
    for molecule in molecules:
        molecule.oneEqvAdded = False
    while True:
        if len(molecules) > MAXLEN:
            #If the reaction gets too crazy, kill.
            raise ReactionTooCrazyError
        places = [(molecule, findPlace(molecule)) for molecule in molecules]
        if not (False in [item[1]==None for item in places]):
            break
        for molecule, place in places:
            if place != None:
                molecules.remove(molecule)
                x = reactAtPlace(molecule, place)
                    
                if not isinstance(x, list):
                    x = [x]
                molecules += x
                if debug:
                    print molecules
    
    if debug:
        print "Results: "
        print smiles(molecules)
    return removeDuplicates(molecules)



def reactWithoutRemoveDuplicates(molecules, findPlace, reactAtPlace):
    if not isinstance(molecules, list):
        return react([molecules], findPlace, reactAtPlace)
    while True:
        places = [(molecule, findPlace(molecule)) for molecule in molecules]
        if not (False in [item[1]==None for item in places]):
            break
        for molecule, place in places:
            if place != None:
                molecules.remove(molecule)
                x = reactAtPlace(molecule, place)
                    
                if not isinstance(x, list):
                    x = [x]
                molecules += x
                if debug:
                    print molecules
    
    if debug:
        print "Results: "
        print smiles(molecules)
    return molecules



"""
Tautomerize
Stuff which looks like ...-C(=C...)-O should actually look like ...-C(-C...)=O
"""
#Does NOT make a duplicate of the list of molecules it takes in.
def tautomerize(moleculeList):

    def reactAtPlace(molecule, place):
        if place==None:
            return molecule
        else:
            target1 = place[0]
            target2 = place[1]
            oxygen = place[2]

            molecule, (target1, target2, oxygen) = listClone(molecule, [target1, target2, oxygen])
            
            #Remove C=C CT-stereocheistry
            target1.eliminateCT()
            target2.eliminateCT()
            target1.eliminateChiral()
            target2.eliminateChiral()
            #Change C-C to single bond
            molecule.changeBond(target1, target2, 1)
            
            #Change C-O to double bond
            molecule.changeBond(target1, oxygen, 2)
            
            return [molecule]


    #Returns a (centralcarbon, alkenecarbon, oxygen) in a tuple. 
    def findPlace(molecule):
        if molecule == None:
            return None
        for atom in molecule.atoms:
            if not (atom.element == 'C'):
                continue
            isAlkene = False
            alkeneCarbon = None
            #check if is alkene
            for neighbor in atom.neighbors:
                if neighbor.element == 'C' and atom.neighbors[neighbor] == 2:
                    isAlkene = True
                    alkeneCarbon = neighbor
            if not isAlkene:
                continue
            #Sanity check this later
            for neighbor in atom.neighbors:
                if neighbor.element == 'O' and atom.neighbors[neighbor] == 1 and len(list(neighbor.neighbors))==1:
                    return (atom, alkeneCarbon, neighbor)
        return None
    
    #if things are crashing, uncomment this line
    #return moleculeList
    
    return reactWithoutRemoveDuplicates(moleculeList, findPlace, reactAtPlace)

     

    #Returns a (centralcarbon, alkenecarbon, oxygen) in a tuple. 
    def findPlace(molecule):
        if molecule == None:
            return None
        for atom in molecule.atoms:
            if not (atom.element == 'C'):
                continue
            isAlkene = False
            alkeneCarbon = None
            #check if is alkene
            for neighbor in atom.neighbors:
                if neighbor.element == 'C' and atom.neighbors[neighbor] == 2:
                    isAlkene = True
                    alkeneCarbon = neighbor
            if not isAlkene:
                continue
            #Sanity check this later
            for neighbor in atom.neighbors:
                if neighbor.element == 'O' and atom.neighbors[neighbor] == 1 and len(list(neighbor.neighbors))==1:
                    return (atom, alkeneCarbon, neighbor)
        return None
    
    #if things are crashing, uncomment this line
    #return moleculeList
    
    return reactWithoutRemoveDuplicates(moleculeList, findPlace, reactAtPlace)



"""Hydrogenation
Candidate reactants: alkenes, alkynes
H2 cat Pd|C in EtOH
Syn addition of an H to each atom in the alkene or alkyne. Go all the way to single bond."""
def hydrogenate(molecules):

    if debug:
        print "Hydrogenating: "
        print smiles(molecules)

    def reactAtPlace(molecule, place):
        if place[0].neighbors[place[1]] == 2:
            #Alkene
            return synAdd(molecule, place[0], place[1], None, None)
        else:
            #Alkyne
            return allTripleAdd(molecule, place[0], place[1], None, None)

    return react(molecules, findAlkeneAndAlkyne, reactAtPlace)



    


"""Hydrohalogenation
HX in CH2Cl2
Candidate reactants: alkenes, alkynes
Adds the X to the Markovnikov-most carbon, and the H to the other carbon. ***Neither syn nor anti*** (because carbocation intermediate).
**UNLESS the alkene is next to a carbonyl (i.e. Michael acceptor). In this case, the halogen is added anti-Markovnikov. We might want to allow users to use this 
reaction, but not allow it to be used when we generate random synthesis problems.
If reacting an alkyne:
if 1eqv specified --> add once
if 2eqv or if excess specified --> add twice
if no quantity specified --> don't let it be a valid reaction? Some sort of feedback to make user specify _how much_ when reacting with alkynes (which is a good habit 
to have) would be nice."""
#halogen is a string
def hydrohalogenate(molecules, halogen):

    if debug:
        print "Hydrohalogenating"
        print smiles(molecules)
        print halogen

    def reactAtPlace(molecule, place): #returns a list of molecules post-reaction at place
        newMolecules = []
        atomicHalogen = Atom(halogen)
        mkvCarbons = markovnikov(place[0], place[1])
        for pairing in mkvCarbons:
            if place[0].neighbors[place[1]] == 2:
                #Double bond
                newMolecules += allAdd(molecule, pairing[0], pairing[1], Atom(halogen), None)
            else:
                #Triple bond
                newMolecules += allTripleAdd(molecule, pairing[0], pairing[1], Atom(halogen), None)

        return newMolecules
    return react(molecules, findAlkeneAndAlkyne, reactAtPlace)


def hydrohalogenate1eq(molecules, halogen):
    
    if debug:
        print "1-eqv. Hydrohalogenating"
        print smiles(molecules)
        print halogen
        
    def findPlace(molecule):
        if hasattr(molecule, 'oneEqvAdded'):
            if molecule.oneEqvAdded:
                return None
        else:
            return None
        if len(findAlkynes(molecule) + findAlkenes(molecule)) != 1:
            return None 
        a = findAlkyne(molecule) 
        if a != None:
            return (a, 3)
        a = findAlkene(molecule)
        if a != None:
            return (a, 2)
        return None
        
    def reactAtPlace(molecule, placeTuple):
        #Check if you are reacting at an alkene or an alkyne.
        if placeTuple == None:
            return []
        place, case = placeTuple
        
        newMolecules = []
        mkvCarbons = markovnikov(place[0], place[1])
        
        for pairing in mkvCarbons:
            if case == 3:
                molecule.oneEqvAdded = True
                newMolecules +=  tripleAdd(molecule, pairing[0], pairing[1], Atom(halogen), None, "cis") + tripleAdd(molecule, pairing[0], pairing[1], Atom(halogen), None, "trans")
            if case == 2:
                molecule.oneEqvAdded = True
                newMolecules +=  allAdd(molecule, pairing[0], pairing[1], Atom(halogen), None)
            
        return newMolecules
        
    return react(molecules, findPlace, reactAtPlace)
    
    
    
    
"""Halogenation
Candidate reactants: alkenes, alkynes
X2 in CH2Cl2
Anti addition of an X to each atom in the alkene.
if 1eqv specified --> add once
if 2eqv or if excess specified --> add twice
if no quantity specified --> don't let it be a valid reaction? Some sort of feedback to make user specify _how much_ when reacting with alkynes (which is a good habit to have) would be nice."""


#halogen is a string
def halogenate(molecules, halogen):
    
    if debug:
        print "Halogenating"
        print smiles(molecules)
        print halogen

    def reactAtPlace(molecule, place):
        atomicHalogen = Atom(halogen)
        atomicHalogen2 = Atom(halogen)
        if place[0].neighbors[place[1]] == 2:
            #Double bond.
            return antiAdd(molecule, place[0], place[1], atomicHalogen, atomicHalogen2)
        else:
            return allTripleAdd(molecule, place[0], place[1], atomicHalogen, atomicHalogen2)
    return react(molecules, findAlkeneAndAlkyne, reactAtPlace)

def halogenate1eq(molecules, halogen):
    
    if debug:
        print "1-eqv. Halogenating"
        print smiles(molecules)
        print halogen
        
    #Reacts one equivalent of X2 with a molecule containing one alkyne and no alkenes.  Will
    #not do anything (e.g. will return the input molecules) if there are multiple alk*nes.
    #Creates a trans alkene.
    #CHANGED: Made this compatible with both alkynes and alkenes.
    def findPlace(molecule):
        if hasattr(molecule, 'oneEqvAdded'):
            if molecule.oneEqvAdded:
                return None
        else:
            return None
        if len(findAlkynes(molecule) + findAlkenes(molecule)) != 1:
            return None 
        a = findAlkyne(molecule) 
        if a != None:
            return (a, 3)
        a = findAlkene(molecule)
        if a != None:
            return (a, 2)
        return None
        
    def reactAtPlace(molecule, placeTuple):
        #Check if you are reacting at an alkene or an alkyne.
        if placeTuple == None:
            return []
        place, case = placeTuple
        if case == 3:
            molecule.oneEqvAdded = True
            return tripleAdd(molecule, place[0], place[1], Atom(halogen), Atom(halogen), 'trans')
        if case == 2:
            molecule.oneEqvAdded = True
            return antiAdd(molecule, place[0], place[1], Atom(halogen), Atom(halogen))
        return []
        
    return react(molecules, findPlace, reactAtPlace)

"""Free-radical hydrohalogenation
Candidate reactants: alkenes
HBr cat ROOR, hv or heat
Adds the X to the anti-Markovnikov-most carbon in the alkene, and the H to the other one. Neither syn nor anti."""

#halogen is a string
def radicalhydrohalogenate(molecules, halogen):
    
    if debug:
        print "Free-radical hydrohalogenating"
        print smiles(molecules)
        print halogen
        
    def findPlace(molecule): #returns one place at which the molecule can react -- e.g. a tuple of atoms, for alkenes/alkynes
        return findAlkene(molecule)
    def reactAtPlace(molecule, place): #returns a list of molecules post-reaction at place
        newMolecules = []
        mkvCarbons = markovnikov(place[0], place[1])
        for pairing in mkvCarbons:
                newMolecules += allAdd(molecule, pairing[0], pairing[1], None, Atom(halogen))
        return newMolecules
    return react(molecules, findPlace, reactAtPlace)



"""
Epoxidation
Candidate reactants: alkenes
mCPBA or PhCO3H or RCO3H, in CH2Cl2
Converts alkene bond to an epoxide. Two possible stereochemical outcomes (up-epoxide or down-epoxide)
"""
def epoxidate(molecules):
    
    if debug:
        print "Epoxidating:"
        print smiles(molecules)

    def epoxAdd(molecule, target1, target2, add1, add2):
        addtarget1 = None
        addtarget2 = None
        antiAdd = False
        #Also does anti-addition, if antiAdd is set to true.
        (molecule, target1, target2, add1, add2, addtarget1, addtarget2) =\
                   duplicateInputs(molecule, target1, target2, add1, add2, addtarget1,
                                   addtarget2)
        (Xmolecule, Xtarget1, Xtarget2, Xadd1, Xadd2, Xaddtarget1, Xaddtarget2) =\
                   duplicateInputs(molecule, target1, target2, add1, add2, addtarget1,
                                   addtarget2)
        add1 = add2
        Xadd1 = Xadd2

        #Set bond orders to single.
        target1.neighbors[target2] = 1
        target2.neighbors[target1] = 1
        Xtarget1.neighbors[Xtarget2] = 1
        Xtarget2.neighbors[Xtarget1] = 1

        bigListOfStuff =\
        ((molecule, add1, target1, addtarget1, target2, target1.CTb, target1.CTa),
        (molecule, add2, target2, addtarget2, target1, target2.CTa, target2.CTb),
        (Xmolecule, Xadd1, Xtarget1, Xaddtarget1, Xtarget2, Xtarget1.CTa, Xtarget1.CTb),
        (Xmolecule, Xadd2, Xtarget2, Xaddtarget2, Xtarget1, Xtarget2.CTb, Xtarget2.CTa))

        for thismolecule, thisAdd, thisTarget, thisAddTarget, otherTarget, ct1, ct2\
                in bigListOfStuff:
            try:
                thismolecule.addAtom(thisAdd, thisTarget, 1)
            except:
                thismolecule.addBond(thisAdd, thisTarget, 1)
            if ct1 != None or ct2 != None:
                thisTarget.newChiralCenter(otherTarget,
                        (thisAdd, ct1, ct2))
            thisTarget.eliminateCT()
        if moleculeCompare(molecule, Xmolecule):
            return [molecule]
        else:
            return [molecule, Xmolecule]
    def findPlace(molecule): #returns one place at which the molecule can react -- e.g. a tuple of atoms, for alkenes/alkynes
        return findAlkene(molecule)
    def reactAtPlace(molecule, place): #returns a list of molecules post-reaction at place
        oxygen = Atom("O")
        return epoxAdd(molecule, place[0], place[1], oxygen, oxygen)
    return react(molecules, findPlace, reactAtPlace)
    
'''
E2 Elimination by tert-butoxide KOC(CH3)3
Removes halide (and a hydrogen on an adjacent carbon) to make a C=C double bond.
If multiple hydrogens can be removed, use Zaitsev's Rule to pick one:
    1) Pick hydrogen that results in more substituted C=C bond.
    2) Bulky, carbon groups trans.
In rings, hydrogen must be anti-periplanar to halide.
'''
def tertButoxide(molecules):
    halogens = ['F', 'Cl', 'Br', 'I']
    if debug:
        print "Applying tertButoxide"

    def findPlace(molecule):
        ans = [] #Used only for complete mode
        #Returns a carbon with a halogen, followed by a carbon with the most suitable H.
        for carbon1 in molecule.atoms:
            if carbon1.element != 'C':
                continue
            #Make sure carbon1 isn't part of a double bond or ketone.
            OK = True
            for neighbor, bo in carbon1.neighbors.items():
                if bo > 1:
                    OK = False
            if not OK:
                continue 
            localHalogens = []
            for neighbor in carbon1.neighbors:
                if (neighbor.element in halogens):
                    localHalogens.append(neighbor)
            if localHalogens == []:
                continue
            #Has to be a carbon, no double or triple bonds, at most 3 things bonded.
            #We don't want to make allenes; and carbonyls have funky behavior that we haven't
            #considered yet.
            neighborCs = [neighbor for neighbor in carbon1.neighbors.keys() if (neighbor.element == 'C'
                and 2 not in neighbor.neighbors.values() and 3 not in neighbor.neighbors.values() and
                len(neighbor.neighbors) <= 3)]
            if len(neighborCs) == 0:
                continue
            #Screw it, let's just try all of them.  Whee, itertools!
            ans += [itertools.product([carbon1], neighborCs, localHalogens)]
        
        return ans
    
    def reactAtPlace(molecule, bigListOfPlaces):
        candidates = []
        for things in bigListOfPlaces:
            #Clone inputs, get on with it.
            Xmolecule, (ClCarbon, HCarbon, Cl) = listClone(molecule, things)
            #Test for epoxides.  We don't deal with epoxides for now.  In reality, attacking an epoxide
            #with KO-tBu results in addition and creation of an ether.
            stop = False
            for neighbor in ClCarbon.neighbors:
                if neighbor in HCarbon.neighbors and neighbor.element == 'O':
                    stop = True
            if stop:
                continue
            #Test for chirality.
            if hasattr(HCarbon, "chiralA") and hasattr(ClCarbon, "chiralA"):
                #Chiral.  We need to consider anti-periplanar.
                #Looking down from the Cl to the other carbon,
                ClsubA, ClsubB = ClCarbon.chiralRingList(Cl, HCarbon)
                #Looking up from the H to the first carbon,
                HsubB, HsubA = HCarbon.chiralRingList(None, ClCarbon)
                #Rings?
                ringList = isInRing(ClCarbon, HCarbon)
                if ringList != None:
                    #The ring carbons must be cis.
                    #ringList: [HCarbon, ..., ClCarbon]
                    if (ClsubA == ringList[-2] and HsubB == ringList[1]) or\
                       (ClsubB == ringList[-2] and HsubA == ringList[1]):
                        pass
                    else:
                        continue
                #Remove chirality, remove XCl, change bond order
                ClCarbon.eliminateChiral()
                HCarbon.eliminateChiral()
                Xmolecule.removeAtom(Cl)
                Xmolecule.changeBond(ClCarbon, HCarbon, 2)
                #Add CTstereo
                try:
                    ClCarbon.newCTCenter(HCarbon, ClsubA, ClsubB)
                    HCarbon.newCTCenter(ClCarbon, HsubA, HsubB)
                except AlleneError:
                    print "Allene error in chiral"
                    continue
                candidates.append((Xmolecule, HCarbon, ClCarbon))
            else:
                #No chirality.
                #May as well set up the double bond now.
                ClCarbon.eliminateChiral()
                HCarbon.eliminateChiral()
                Xmolecule.removeAtom(Cl)
                Xmolecule.changeBond(ClCarbon, HCarbon, 2)
                #Rings?
                ringList = isInRing(ClCarbon, HCarbon)
                if ringList != None:
                    #Make rings cis.
                    ClsubA = ringList[-2]
                    ClsubB = None
                    for neighbor in ClCarbon.neighbors:
                        if neighbor != ClsubA and neighbor != HCarbon:
                            ClsubB = neighbor
                    HsubB = ringList[1]
                    HsubA = None
                    for neighbor in HCarbon.neighbors:
                        if neighbor != HsubB and neighbor != ClCarbon:
                            HsubA = neighbor
                    try:
                        ClCarbon.newCTCenter(HCarbon, ClsubA, ClsubB)
                        HCarbon.newCTCenter(ClCarbon, HsubA, HsubB)
                    except AlleneError:
                        print "Allene error in ring, no chiral"
                        continue
                    candidates.append((Xmolecule, HCarbon, ClCarbon))
                else:
                    #No rings.  Make both cases.
                    Clsubs = []
                    for neighbor in ClCarbon.neighbors:
                        if neighbor != HCarbon:
                            Clsubs.append(neighbor)
                    while len(Clsubs) < 2:
                        Clsubs.append(None)
                    Hsubs = []
                    for neighbor in HCarbon.neighbors:
                        if neighbor != ClCarbon:
                            Hsubs.append(neighbor)
                    while len(Hsubs) < 2:
                        Hsubs.append(None)
                    Clsubs2 = [[],[]]
                    Hsubs2 = [[],[]]
                    #Does each carbon have exactly one other substituent?  If so, cis/trans stereochem
                    #becomes important.
                    if len(ClCarbon.neighbors) == 2 and len(HCarbon.neighbors) == 2:
                        #Return only the trans molecule.
                        try:
                            ClCarbon.newCTCenter(HCarbon, Clsubs[0], Clsubs[1])
                            HCarbon.newCTCenter(ClCarbon, Hsubs[0], Hsubs[1])
                        except AlleneError:
                            print "Allene error in noring"
                            continue
                        candidates.append((Xmolecule, HCarbon, ClCarbon))
                    else:
                        Xmolecule2, (ClCarbon2, HCarbon2, Clsubs2[0], Clsubs2[1], Hsubs2[0], Hsubs2[1]) =\
                            listClone(Xmolecule, (ClCarbon, HCarbon, Clsubs[0], Clsubs[1], Hsubs[0], Hsubs[1]))
                        try:
                            ClCarbon.newCTCenter(HCarbon, Clsubs[0], Clsubs[1])
                            HCarbon.newCTCenter(ClCarbon, Hsubs[0], Hsubs[1])
                            candidates.append((Xmolecule, HCarbon, ClCarbon))
                        except AlleneError:
                            print "Allene error in noringc/t"
                        try:
                            ClCarbon2.newCTCenter(HCarbon2, Clsubs2[1], Clsubs2[0])
                            HCarbon2.newCTCenter(ClCarbon2, Hsubs2[0], Hsubs2[1])
                            candidates.append((Xmolecule2, HCarbon2, ClCarbon2))
                        except AlleneError:
                            print "Allene error in noringc/t"

        #Now, prune the products
        #Find the most substituted products, and keep only those.
        maxSub = 0
        for Xmolecule, c1, c2 in candidates:
            maxSub = max(maxSub, len(c1.neighbors)+len(c2.neighbors))
        maxSubCandidates = []
        for Xmolecule, c1, c2 in candidates:
            if len(c1.neighbors)+len(c2.neighbors) == maxSub:
                maxSubCandidates.append(Xmolecule)
        if debug:
            print maxSubCandidates
        return maxSubCandidates

        
    def complete(molecules):
        #Complete testing.
        if len(molecules) > MAXLEN:
            raise ReactionTooCrazyError
        out = []
        #Note to self: do not modify molecules.  You need it for returning at the end.
        for molecule in molecules:
            places = findPlace(molecule)
            if places == []:
                out += [molecule]
                continue
            thisOut = []
            for place in places:
                ans = reactAtPlace(molecule, place)
                if debug:
                    for thing in ans:
                        verify(thing)
                thisOut+=ans
            if thisOut == []:
                out += [molecule]
            else:
                out += thisOut
        if out==molecules:
            return molecules
        else:
            if debug:
                print out
                print molecules
        return complete(removeDuplicates(out))
    if debug:
        a = complete(molecules)
        print "Result of tBut: " +str(smiles(a))
        return a
    return removeDuplicates(complete(molecules))
    
    
def listClone(molecule, atomList):
    #Returns properly-connected deepcopies of molecule and list of atoms.
    Xmolecule = copy.deepcopy(molecule)
    newAtomList = []
    for atom in atomList:
        if atom == None:
            newAtomList.append(None)
        else:
            newAtomList.append(Xmolecule.atoms[molecule.atoms.index(atom)])
    return Xmolecule, newAtomList
          


"""
Hydration
Candidate reactants: alkenes
H2SO4 (or other acid) in ROH, where R can also be H
If alkene: Adds an OR to the Markovnikov carbon of the alkene, and an H to the anti-Markovnikov carbon. Neither syn nor anti, since carbocation.
If alkyne and H2O: Form a ketone or aldehyde, placing the O at the Markovnikov carbon.
If alkyne and ROH: form an enolate-ether, I think...     or form two ethers.

For alkynes to react, usually also mention "HgSO4 accels."
"""

#When there is an other-molecule:
    #If any incoming molecules can react with themselves, use that as the product.
    #Only if none of the molecules can react with themselves, react them against each other in all possible ways.
        #Afterwards, check the resulting molecule for self-reactivity.

def acidhydrate(molecules, others, alkynesOk = False):


    def findPlaces1(molecule):
        if alkynesOk:
            return findAlkynes(molecule)
        else:
            return findAlkenes(molecule)
    def findPlaces2(molecule):
        return findHydroxyls(molecule)

    #Place 1 is an alkene or an alkyne (tuple of atoms)
    #Place 2 is an oxygen connected to 1 or 0 neighbors
    def reactAtPlaces(molecule1, molecule2, place1, place2):
        
        #Use addMolecule(self, molecule, foreignTarget, selfTarget, bo)
        if place1[0].neighbors[place1[1]] == 2: #if is alkene:
            if debug:
                print "Case 1: alkene in acidhydrate"
            newMolecules = []
            mkvCarbons = markovnikov(place1[0], place1[1])
            for pairing in mkvCarbons:
                newMolecules += allAdd(molecule1, pairing[0], pairing[1], molecule2, None, place2, None)
            return newMolecules

        elif place1[0].neighbors[place1[1]] == 3: #if is alkyne:
            if debug:
                print "Case 2: alkyne in acidhydrate"
            
            if len(list(place2.neighbors)) == 0: #if adding water:
                if debug:
                    print "Case 3: ...with water"
                #Going to need to write a custom function here, borrowed from allTripleAdd.
                #Make the alkyne bond a single bond
                #Add a double-bond-O to the Markovnikov carbon
                #If each carbon is equally Markovnikov, do some duplication-hacks and add to both
                #Make a double bond between each Markovnikov carbon and the O of place2
                
                newMolecules = []
                mkvCarbons = markovnikov(place1[0], place1[1])
                for pairing in mkvCarbons:
                    newMolecules += carbonylAdd(molecule1, pairing[0], pairing[1])
                return newMolecules

                
            elif len(list(place2.neighbors)) == 1: #if is alcohol:
                if debug:
                    print "Case 4: ...with alcohol"
                #Not sure if chemically correct - FS.  Definitely not supposed to be
                #allTripleAdd - that results in an alkane.
                
                newMolecules = []
                mkvCarbons = markovnikov(place1[0], place1[1])
                for pairing in mkvCarbons:
                    newMolecules += tripleAdd(molecule1, pairing[0], pairing[1], molecule2, None, 'trans', place2, None, )
                return newMolecules
        else:
            print "Error: findAlkenes, findAlkynes returning non-alkene and/or non-alkyne"
            raise StandardError

        
        
        
    return twoReact(copy.deepcopy(molecules), copy.deepcopy(others), findPlaces1, findPlaces2, reactAtPlaces)


#NOTE: findPlaces methods passed into this method MUST return lists
def twoReact(molecules, others, findPlaces1, findPlaces2, reactAtPlaces):
    if not isinstance(molecules, list):
        return twoReact([molecules], others, findPlaces1, findPlaces2, reactAtPlaces)
    if not isinstance(others, list):
        return twoReact(molecules, [others], findPlaces1, findPlaces2, reactAtPlaces)
    output = []
    molecules1 = [] #molecules which are capable of playing role 1
    molecules2 = [] #molecules which are capable of playing role 2
    for molecule in molecules+others:
        candidates1 = [x for x in findPlaces1(molecule) if x != None] #places in molecule which can react as role 1
        candidates2 = [x for x in findPlaces2(molecule) if x != None] #places in molecule which can react as role 2
        if len(candidates1) != 0:
            if len(candidates2) != 0:
                #self-react and add to list
                output += [item for sublist in [
                    reactAtPlaces(molecule, molecule, locus1, locus2)
                        for locus1 in findPlaces1(molecule)
                        for locus2 in findPlaces2(molecule)
                    ]
                    for item in sublist]

        if len(candidates1) != 0:
            molecules1 += [molecule]

        if len(candidates2) != 0:
            molecules2 += [molecule]

    #If this is true, then no molecule reacted with itself. You may proceed to reacting the molecules in molecules1 and molecules2 with each other.
    if len(output) == 0:
        output = [item for sublist in [reactAtPlaces(molecule1, molecule2, locus1, locus2) for molecule1 in molecules1 for molecule2 in molecules2 for locus1 in findPlaces1(molecule1) for locus2 in findPlaces2(molecule2)] for item in sublist]

    if debug:
        print "Results: "
        print smiles(output)
        
    return removeDuplicates(output)





"""
Halohydration
Candidate reactants: alkenes
X2 in ROH, where R can also be H
Adds an OR to the Markovnikov carbon of the alkene, and an X to the anti-Markovnikov carbon. Anti.
"""
def halohydrate(molecules, others, halogen):
    
    if debug:
        print "Halohydrating"
        print smiles(molecules)
        print smiles(others)
        print halogen

    atomicHalogen = Atom(halogen)
    def findPlaces1(molecule):
        return findAlkenes(molecule)
    def findPlaces2(molecule):
        return findHydroxyls(molecule)

    #Place 1 is an alkene (tuple of atoms)
    #Place 2 is an oxygen connected to 1 or 0 neighbors
    def reactAtPlaces(molecule1, molecule2, place1, place2):
        newMolecules = []
        mkvCarbons = markovnikov(place1[0], place1[1])
        for pairing in mkvCarbons:
            newMolecules += antiAdd(molecule1, pairing[0], pairing[1], molecule2, atomicHalogen, place2, None)
        return newMolecules
        
    return twoReact(copy.deepcopy(molecules), copy.deepcopy(others), findPlaces1, findPlaces2, reactAtPlaces)



"""
Hydroboration
Candidate reactants: alkenes, alkynes
BH3 in THF, then NaOH and H2O2
If alkene: Adds an OH to the anti-Markovnikov carbon, then adds an H to the other one. Syn addition.
If alkyne: Form a ketone or aldehyde, placing the O at the anti-Markovnikov carbon.
"""
#Both steps at once
def hydroborate(molecules):
    
    if debug:
        print "Hydroborating"
        print smiles(molecules)
        
        
    def reactAtPlace(molecule, place): #returns a list of molecules post-reaction at place
        newMolecules = []
        oxy = Atom("O")
        mkvCarbons = markovnikov(place[0], place[1])
        if place[0].neighbors[place[1]] == 2: #Alkene
            for pairing in mkvCarbons:
                newMolecules += synAdd(molecule, pairing[0], pairing[1], None, oxy)
        else: #Alkyne
            for pairing in mkvCarbons:
                newMolecules += carbonylAdd(molecule, pairing[1], pairing[0])
        return newMolecules
    
    return react(molecules, findAlkeneAndAlkyne, reactAtPlace)

#Only BH3 in THF
def hydroborate1(molecules):
    
    if debug:
        print "Hydroborating 1"
        print smiles(molecules)
        
    def findOkAlkene(molecule):
        #any alkene which already has a borane attached is invalid
        for atom in molecule.atoms:
            if not (atom.element == 'C'):
                continue
            r = None
            for neighbor in atom.neighbors:
                if neighbor.element == 'C' and atom.neighbors[neighbor] == 2:
                    r = (atom, neighbor)
            if r != None:
                if not ("B" in [n.element for n in r[0].neighbors] or "B" in [n.element for n in r[1].neighbors]):
                    return r
        return None
        
    def findPlace(molecule):
        x = findOkAlkene(molecule)
        if x == None:
            return findAlkyne(molecule)
        else:
            return x
        
    
    def reactAtPlace(molecule, place): #returns a list of molecules post-reaction at place
        newMolecules = []
        boron = Atom("B")
        mkvCarbons = markovnikov(place[0], place[1])
        if place[0].neighbors[place[1]] == 2: #Alkene
            for pairing in mkvCarbons:
                newMolecules += synAdd(molecule, pairing[0], pairing[1], None, boron)
        else: #Alkyne
        #NO. NOT CARBONYL-ADD.
            for pairing in mkvCarbons:
                newMolecules += tripleAdd(molecule, pairing[1], pairing[0], boron, None, "cis")
        return newMolecules
    
    return react(molecules, findPlace, reactAtPlace)

#Subsequent NaOH and H2O2 step
def hydroborate2(molecules):
    
    if debug:
        print "Hydroborating 2"
        print smiles(molecules)
        
    def BtoO(molecule):
        for atom in molecule.atoms:
            if atom.element == "B" and len(list(atom.neighbors)) < 2:
                atom.element = "O"
        return molecule
    newMolecules = copy.deepcopy(molecules)
    return removeDuplicates([BtoO(molecule) for molecule in newMolecules])


"""
Dihydroxylation (Upjohn dihydroxylation)
Candidate reactants: alkenes
cat. OsO4 in NMO and acetone or H2O
Syn addition of two OH groups to each carbon.
"""
def dihydroxylate(molecules):
    
    if debug:
        print "Dihydroxylating"
        print smiles(molecules)
        
    def reactAtPlace(molecule, place):
        return synAdd(molecule, place[0], place[1], Atom("O"), Atom("O"))
    return react(molecules, findAlkene, reactAtPlace)



"""
Ozonolysis
Candidate reactants: alkenes
O3 in CH2Cl2, with Me2S or Zn
Adds two oxygens, splitting alkene bond, producing carbonyls.
"""
def ozonolyse(molecules):
    
    if debug:
        print "Ozonolysing"
        print smiles(molecules)
        
        
    def reactAtPlace(molecule, place):
        #Break the double bond
        molecule.changeBond(place[0], place[1], 0)
        #Add two double bonds to two new oxygens
        molecule.addAtom(Atom("O"), place[0], 2)
        molecule.addAtom(Atom("O"), place[1], 2)
        #Destroy CT stereochemistry
        place[0].eliminateCT()
        place[1].eliminateCT()
        place[0].eliminateChiral()
        place[1].eliminateChiral()
        #Splice the molecule
        return splice(molecule)
    return react(copy.deepcopy(molecules), findAlkene, reactAtPlace)




"""
Lindlar reduction
Candidate reactants: alkynes
H2, cat. Lindlar
Produces the cis alkene from an alkyne. Adds two Hs.
"""
def lindlar(molecules):
    
    if debug:
        print "Lindlar-ing"
        print smiles(molecules)
        
    def reactAtPlace(molecule, place):
        return tripleAdd(molecule, place[0], place[1], None, None, 'cis')
    return react(molecules, findAlkyne, reactAtPlace)



"""
Sodium-ammonia reduction
Candidate reactants: alkynes
Na, in NH3 (L)
Produces the trans alkene from an alkyne. Adds two Hs.
"""
def sodiumAmmonia(molecules):
    
    if debug:
        print "Sodium-ammonia-reducing"
        print smiles(molecules)
        
    def reactAtPlace(molecule, place):
        return tripleAdd(molecule, place[0], place[1], None, None, 'trans')
    return react(molecules, findAlkyne, reactAtPlace)




"""
Alkyne deprotonation to acetylide
Candidate reactants: alkynes, in which one end is an H
NaNH2 in NH3
Produces an acetylide ion. Removes the H+, resulting in a negative charge.
"""
def alkyneDeprotonate(molecules):
    
    if debug:
        print "Alkyne-deprotonating"
        print smiles(molecules)
        
    def findPlace(molecule):
        if debug:
            print molecule
        return findHydrogenAlkyne(molecule)
    def reactAtPlace(molecule, place):
        if place == None:
            return []
        
        place.charge = -1
        return [molecule]
        
    return react(copy.deepcopy(molecules), findPlace, reactAtPlace)





"""
Reactants: one acetylide ion, one molecule of the form R-CH2-Br or R-CH2-I
No additional reagents needed (can drag one molecule onto the other?)
Get rid of the X in R-CH2-X.
Attach the bare negative end of the acetylide (R-C#Cminus) to the R-CH2-.
The result should look like R-CH2-C#C-R".
"""
HALOGENS = ["Br","I", "Cl"]
def acetylideAdd(molecules, others):
    
    if debug:
        print "Acetylide-adding"
        print smiles(molecules)
        print smiles(others)
    
    def findPlaces1(molecule):
        #findAlkyneCarbanions(molecule)
        places = []
        for atom in molecule.atoms:
            if atom.charge == -1 and atom.element == "C" and (True in [other.neighbors[atom]==3 for other in list(atom.neighbors) if other.element == "C"]):
                places += [atom]
        return places
    def findPlaces2(molecule):
        #findHalogenCarbons(molecule)
        places = []
        countedOxygens = []
        for atom in molecule.atoms:
            #Primary alkyl halides
            if atom.element == "C" and len(list(atom.neighbors))==2 and (True in [(other.element in HALOGENS) for other in list(atom.neighbors)]):
                places += [atom]
            #Epoxides
            if atom.element == 'C':
                for neighbor in atom.neighbors:
                    if neighbor.element == 'O' and neighbor not in countedOxygens:
                        for otherN in atom.neighbors:
                            if otherN in neighbor.neighbors.keys():
                                #carbon, other carbon, oxygen
                                places += [(atom, otherN, neighbor)]
                                countedOxygens.append(neighbor)
                
        return places

    #Place 1 is a negatively charged carbon atom
    #Place 2 is a carbon connected to at least one halogen, or an epoxide tuple
    def reactAtPlaces(molecule1, molecule2, place1, place2):
        (molecule1, place1, unused0, unused1, unused2, unused3, unused4)=\
               duplicateInputs(molecule1, place1, place1, None, None, None, None)
        if isinstance(place2, tuple):
            #Epoxide.  I just realized that this can also be implemented with antiAdd, but oh well.
            #Do some object acrobatics to make proper deepcopies.
            oxIndex = molecule2.atoms.index(place2[2])
            (molecule2, carbon, otherCarbon, unused0, unused1, unused3, unused2) =\
                duplicateInputs(molecule2, place2[0], place2[1], None, None, None, None)
            oxygen = molecule2.atoms[oxIndex]
            mkvPairs = markovnikov(carbon, otherCarbon)
            if len(mkvPairs) == 2:
                (Xmolecule1, Xplace1, unused0, unused1, unused2, unused3, unused4)=\
                    duplicateInputs(molecule1, place1, place1, None, None, None, None)
                (Xmolecule2, addingC, alcoholC, unused0, unused1, unused3, unused2) =\
                    duplicateInputs(molecule2, mkvPairs[1][0], mkvPairs[1][1], None, None, None, None)
                Xoxygen = Xmolecule2.atoms[oxIndex]
                #I apoplogize for the rather confusing naming - ran out of good name ideas.
                stuff = ((molecule1, place1, molecule2, mkvPairs[0][0], mkvPairs[0][1], oxygen),
                     (Xmolecule1, Xplace1, Xmolecule2, addingC, alcoholC, Xoxygen))
            else:
                stuff = ((molecule1, place1, molecule2, mkvPairs[0][0], mkvPairs[0][1], oxygen))

            for (aceMol, aceAtom, epxMol, addC, alcC, oxygen) in stuff:
                #Break the epoxide bond.
                epxMol.changeBond(addC, oxygen, 0)
                #Remove epoxide stereochem!
                viewFromOx = addC.chiralCWlist(oxygen)
                addC.eliminateChiral()
                #alcC keeps its stereochemistry
                #Add the acetylene to the more substituted carbon.
                aceAtom.charge = 0
                epxMol.addMolecule(aceMol, aceAtom, addC, 1)
                #Add new stereochemistry - inversion.
                addC.newChiralCenter(aceAtom, (viewFromOx[0], viewFromOx[2], viewFromOx[1]))
            if len(mkvPairs) == 2:
                if moleculeCompare(molecule2, Xmolecule2):
                    return [molecule2]
                else:
                    return [molecule2, Xmolecule2]
            else:
                return [molecule2]
         
        else:
            #Alkyl halide
            #Duplicate inputs
            (molecule1, place1, unused0, unused1, unused2, unused3, unused4)=\
                   duplicateInputs(molecule1, place1, place1, None, None, None, None)
            (molecule2, place2, unused0, unused1, unused2, unused3, unused4)=\
                   duplicateInputs(molecule2, place2, place2, None, None, None, None)
           
            #Remove negative charge
            place1.charge = 0
            #Find halogen
            halogen = None
            for atom in list(place2.neighbors):
                if atom.element in HALOGENS:
                    halogen = atom
            if halogen == None:
                print "Error: your HalogenCarbons method isn't working as intended"
                raise StandardError
            #Remove halogen
            molecule2.removeAtom(halogen)
            #Add a single bond between the two carbons
            molecule1.addMolecule(molecule2, place2, place1, 1)
        
        return [molecule1]
        
        
    return twoReact(copy.deepcopy(molecules), copy.deepcopy(others), findPlaces1, findPlaces2, reactAtPlaces)










#DO NOT DELETE. This is referred to by code in other files.
#Makes C-OH, methanol.
c72 = Atom("C")
o73 = Atom("O")
methanol = Molecule(c72)
methanol.addAtom(o73, c72, 1)

#DO NOT DELETE. This is referred to by code in other files.
#Makes C-C-OH, ethanol.
c64 = Atom("C")
c65 = Atom("C")
o66 = Atom("O")
ethanol = Molecule(c64)
ethanol.addAtom(c65, c64, 1)
ethanol.addAtom(o66, c64, 1)


if __name__ == '__main__':
    #Makes     C-C-C<C
    #          |   |
    #        O-C=C-N

    #Makes     C-C-C>C
    #          |   |
    #        O-C=C-N
    c1 = Atom("C")
    mol = Molecule(c1)
    c2 = Atom("C")
    n1 = Atom("N")
    mol.addAtom(c2, c1, 2)
    mol.addAtom(n1, c2, 1)
    o1 = Atom("O")
    mol.addAtom(o1, c1, 1)

    c3 = Atom("C")
    mol.addAtom(c3, n1, 1)
    c4 = Atom("C")
    c5 = Atom("C")
    mol.addAtom(c4, c3, 1)
    mol.addAtom(c5, c3, 1)
    c6 = Atom("C")
    mol.addAtom(c6, c5, 1)
    mol.addBond(c6, c1, 1)
    c3.newChiralCenter(n1, (c4, None, c5))
    c1.newCTCenter(c2, o1, c6)
    c2.newCTCenter(c1, n1, None)



    #Makes C\   /Cl
    #        C=C
    #     Cl/   \Br
    c10 = Atom("C")
    CTmol = Molecule(c10)
    c11 = Atom("C")
    CTmol.addAtom(c11, c10, 2)
    c12 = Atom("C")
    CTmol.addAtom(c12, c10, 1)
    cL1 = Atom("Cl")
    CTmol.addAtom(cL1, c10, 1)
    cL2 = Atom("Cl")
    CTmol.addAtom(cL2, c11, 1)
    br10 = Atom("Br")
    CTmol.addAtom(br10, c11, 1)
    c10.newCTCenter(c11, cL1, c12)
    c11.newCTCenter(c10, cL2, br10)

    #Makes C\   /Br
    #        C=C
    #     C1/   \Cl
    c15 = Atom("C")
    CTmol2 = Molecule(c15)
    c16 = Atom("C")
    CTmol2.addAtom(c16, c15, 2)
    c17 = Atom("C")
    CTmol2.addAtom(c17, c15, 1)
    cL5 = Atom("Cl")
    CTmol2.addAtom(cL5, c15, 1)
    cL6 = Atom("Cl")
    CTmol2.addAtom(cL6, c16, 1)
    br15 = Atom("Br")
    CTmol2.addAtom(br15, c16, 1)
    c15.newCTCenter(c16, cL5, c17)
    c16.newCTCenter(c15, br15, cL6)

    #Makes  C\ /C-C
    #         C
    #      Br/ \H
    c20 = Atom("C")
    chiralMol1 = Molecule(c20)
    c23 = Atom("C")
    chiralMol1.addAtom(c23, c20, 1)
    br20 = Atom("Br")
    chiralMol1.addAtom(br20, c20, 1)
    c21 = Atom("C")
    chiralMol1.addAtom(c21, c20, 1)
    c22 = Atom("C")
    chiralMol1.addAtom(c22, c21, 1)
    c20.newChiralCenter(c21, (None, br20, c23))

    c30 = Atom("C")
    chiralMol2 = Molecule(c30)
    c33 = Atom("C")
    chiralMol2.addAtom(c33, c30, 1)
    br30 = Atom("Br")
    chiralMol2.addAtom(br30, c30, 1)
    c31 = Atom("C")
    chiralMol2.addAtom(c31, c30, 1)
    c32 = Atom("C")
    chiralMol2.addAtom(c32, c31, 1)
    c30.newChiralCenter(c31, (None, c33, br30))

    #Makes C/C=C
    c40 = Atom("C")
    c41 = Atom("C")
    mol4 = Molecule(c40)
    mol4.addAtom(c41, c40, 2)
    c42 = Atom("C")
    mol4.addAtom(c42, c40, 1)
    c40.newCTCenter(c41, None, c42)
    c41.newCTCenter(c40, None, None)

    #Makes C\C=C
    c45 = Atom("C")
    c46 = Atom("C")
    mol4alt = Molecule(c45)
    mol4alt.addAtom(c46, c45, 2)
    c47 = Atom("C")
    mol4alt.addAtom(c47, c45, 1)
    c45.newCTCenter(c46, c47, None)
    c46.newCTCenter(c45, None, None)

    #        c50
    #Makes C-C<Cl
    #     /   \
    #   C<C   C>Br c51
    #      \C/
    c50 = Atom("C")
    c51 = Atom("C")
    c52 = Atom("C")
    c53 = Atom("C")
    c54 = Atom("C")
    cycPentMol = Molecule(c50)
    cycPentMol.addAtom(c51, c50, 1)
    cycPentMol.addAtom(c52, c51, 1)
    cycPentMol.addAtom(c53, c52, 1)
    cycPentMol.addAtom(c54, c53, 1)
    cycPentMol.addBond(c54, c50, 1)
    cl50 = Atom("Cl")
    cycPentMol.addAtom(cl50, c50, 1)
    c50.newChiralCenter(c54, (cl50, c51, None))
    br50 = Atom("Br")
    cycPentMol.addAtom(br50, c51, 1)
    c51.newChiralCenter(c50, (c52, br50, None))
    c55 = Atom("C")
    cycPentMol.addAtom(c55, c53, 1)
    c53.newChiralCenter(c52, (c55, None, c54))
    
    #Makes C-C>Cl
    #     /   \
    #   C<C   C>Br c51
    #      \C/
    c150 = Atom("C")
    c151 = Atom("C")
    c152 = Atom("C")
    c153 = Atom("C")
    c154 = Atom("C")
    cycPentMol2 = Molecule(c150)
    cycPentMol2.addAtom(c151, c150, 1)
    cycPentMol2.addAtom(c152, c151, 1)
    cycPentMol2.addAtom(c153, c152, 1)
    cycPentMol2.addAtom(c154, c153, 1)
    cycPentMol2.addBond(c154, c150, 1)
    cl150 = Atom("Cl")
    cycPentMol2.addAtom(cl150, c150, 1)
    c150.newChiralCenter(c154, (c151, cl150, None))
    br150 = Atom("Br")
    cycPentMol2.addAtom(br150, c151, 1)
    c151.newChiralCenter(c150, (c152, br150, None))
    c155 = Atom("C")
    cycPentMol2.addAtom(c155, c153, 1)
    c153.newChiralCenter(c152, (c155, None, c154))

    #Makes C-C#C-C
    c60 = Atom("C")
    c61 = Atom("C")
    c62 = Atom("C")
    c63 = Atom("C")
    propyne = Molecule(c60)
    propyne.addAtom(c61, c60, 3)
    propyne.addAtom(c62, c61, 1)
    propyne.addAtom(c63, c60, 1)

    #Makes C#C
    c67 = Atom("C")
    c68 = Atom("C")
    ethylene = Molecule(c67)
    ethylene.addAtom(c68, c67, 3)

    #Makes C-C-Br
    c69 = Atom("C")
    br70 = Atom("Br")
    c71 = Atom("C")
    bromoethane = Molecule(c69)
    bromoethane.addAtom(c71, c69, 1)
    bromoethane.addAtom(br70, c69, 1)

    #       F
    #Makes CCC
    #       CC
    c80 = Atom("C")
    c81 = Atom("C")
    c82 = Atom("C")
    c83 = Atom("C")
    c84 = Atom("C")
    c85 = Atom("C")
    c86 = Atom("C")
    c87 = Atom("C")
    c88 = Atom("C")
    c89 = Atom("C")
    f80 = Atom("F")
    f81 = Atom("F")
    ringTest1 = Molecule(c80) #Tail
    ringTest1.addAtom(c81, c80, 1)
    ringTest1.addAtom(c82, c81, 1)
    ringTest1.addAtom(c83, c82, 1)
    ringTest1.addAtom(c84, c83, 1)
    ringTest1.addBond(c84, c81, 1)
    ringTest1.addAtom(f80, c81, 1)
    c81.newChiralCenter(c80, (c82, c84, f80))

    ringTest2 = Molecule(c85) #Ring
    ringTest2.addAtom(c86, c85, 1) #Tail
    ringTest2.addAtom(c87, c85, 1)
    ringTest2.addAtom(c88, c87, 1)
    ringTest2.addAtom(c89, c88, 1)
    ringTest2.addBond(c89, c85, 1)
    ringTest2.addAtom(f81, c85, 1)
    c85.newChiralCenter(c86, (c89, c87, f81))



    c90 = Atom("C")
    c91 = Atom("C")
    c92 = Atom("C")
    c93 = Atom("C")
    c94 = Atom("C")
    o95 = Atom("O")
    ring = Molecule(c90)
    ring.addAtom(c91, c90, 2)
    ring.addAtom(c92, c91, 1)
    ring.addAtom(c93, c92, 1)
    ring.addAtom(c94, c93, 1)
    ring.addBond(c90, c94, 1)

    ring.addAtom(o95, c94, 1)
    ring.addBond(c90, o95, 1)

    
    c90.newCTCenter(c91, c94, None)
    c91.newCTCenter(c90, c92, None)

    
    print reduceChirality([cycPentMol, cycPentMol2])
