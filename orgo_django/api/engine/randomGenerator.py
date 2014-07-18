from helperFunctions import *
import random

def randomStart(endProb=0.3, maxBranchLength=15,
                alkyneProb=0.2, alkeneProb=0.4,
                BrProb=0.1, ClProb=0.1, OHProb=0.1, BranchProb=0.05,
                forceTerminalAlkyne = False):
    #Want more alkynes?  You need to make fewer substituents.
    for subsProb in (BrProb, ClProb, OHProb, BranchProb):
        subsProb *= 1-2*alkyneProb
    lastAtom = Atom("C")
    frontAtom = lastAtom
    mol = Molecule(lastAtom)
    if forceTerminalAlkyne:
        lastAtom = Atom("C")
        mol.addAtom(lastAtom, frontAtom, 3)
    while (random.random() < 1.0-endProb or len(mol.atoms) < 3) and len(mol.atoms) < maxBranchLength:
        switcher = random.random()
        if switcher < 0.8:
            newMol, thisAtom, nextAtom = randC(BrProb, ClProb, OHProb, BranchProb)
        elif switcher < 0.9:
            newMol, thisAtom, nextAtom = randRing(5, BrProb, ClProb, OHProb)
        else:
            newMol, thisAtom, nextAtom = randRing(6, BrProb, ClProb, OHProb)

        switcher = random.random()
        if switcher < alkyneProb and len(thisAtom.neighbors) == 0\
           and lastAtom.totalBondOrder() == 1:
            #Make a triple bond.
            mol.addMolecule(newMol, thisAtom, lastAtom, 3)
        elif switcher < alkeneProb+alkyneProb and len(thisAtom.neighbors) <= 1\
             and lastAtom.totalBondOrder() <= 2\
             and lastAtom.findAlkeneBond() == None\
             and lastAtom.totalBondOrder() > 0:
            #Make a double bond.  Watch out for cis/trans.
            #Each double bond must have a cis/trans specification.
            mol.addMolecule(newMol, thisAtom, lastAtom, 2)
            otherNeighbors = []
            for neighbor in lastAtom.neighbors:
                if neighbor != thisAtom:
                    otherNeighbors.append(neighbor)
            if len(otherNeighbors) == 1:
                lastAtom.newCTCenter(thisAtom, otherNeighbors[0],
                                     None)
            else:
                lastAtom.newCTCenter(thisAtom, otherNeighbors[0],
                                     otherNeighbors[1])
            #Flag thisAtom so that we add stereo to it, later.
            thisAtom.makeCTFlag = True
        else:
            #Make single bond
            mol.addMolecule(newMol, thisAtom, lastAtom, 1)
        fixStereo(mol, thisAtom, lastAtom)
        lastAtom = thisAtom
        thisAtom = nextAtom
    #Not quite done yet.  We need to fix stereochemistry on the last atom.
    fixStereo(mol, None, lastAtom) #May be a slight problem with chirality?
    return mol, frontAtom, None

def randC(ClProb, BrProb, OHProb, BranchProb):
    #Makes a random carbon atom with some substituents
    c = Atom("C")
    mol = Molecule(c)
    #Add up to 2 substituents
    for i in xrange(2):
        switcher = random.random()
        if switcher < ClProb:
            newS = Atom("Cl")
            mol.addAtom(newS, c, 1)
        elif switcher < ClProb+BrProb:
            newS = Atom("Br")
            mol.addAtom(newS, c, 1)
        elif switcher < ClProb+BrProb+OHProb:
            newS = Atom("O")
            mol.addAtom(newS, c, 1)
        elif switcher < ClProb+BrProb+OHProb+BranchProb:
            print "Branching"
            newS, frontAtom, notused = randomStart()
            mol.addMolecule(newS, frontAtom, c, 1)
    return mol, c, c

def randRing(noCs, BrProb, ClProb, OHProb):
    initAtom = Atom("C")
    mol = Molecule(initAtom)
    oldAtom = initAtom
    for i in xrange(noCs - 1):
        newAtom = Atom('C')
        mol.addAtom(newAtom, oldAtom, 1)
        oldAtom = newAtom
    #Now, close the loop!
    mol.addBond(initAtom, oldAtom, 1)
    outAtom = mol.atoms[random.randint(1,len(mol.atoms)-1)]
    for atom in mol.atoms:
        #Don't interfere with the larger structure.
        if atom == initAtom or atom == outAtom:
            continue
        #Add some random substituents.
        switcher = random.random()
        if switcher < BrProb:
            addAtom = Atom("Br")
        elif switcher < BrProb + ClProb:
            addAtom = Atom("Cl")
        elif switcher < BrProb + ClProb + OHProb:
            addAtom = Atom("O")
        else:
            continue
        if len(atom.neighbors) < 2:
            continue
        #Do a coin flip to determine chirality.
        if random.random() < 0.5:
            atom.newChiralCenter(addAtom, (atom.neighbors.keys()[0], atom.neighbors.keys()[1], None))
        else:
            atom.newChiralCenter(addAtom, (atom.neighbors.keys()[1], atom.neighbors.keys()[0], None))
        mol.addAtom(addAtom, atom, 1)
    
    return mol, initAtom, outAtom

def fixStereo(mol, thisAtom, lastAtom):
    #Look to see if the *last* piece we added requires stereochem
    otherC = lastAtom.findAlkeneBond()
    if hasattr(lastAtom, "makeCTFlag"):
        del(lastAtom.makeCTFlag)
        #Find the last neighbor of lastAtom (other than thisAtom
        #and otherC).  None means Hydrogen.
        otherN = None
        for neighbor in lastAtom.neighbors:
            if neighbor != otherC and neighbor != thisAtom:
                otherN = neighbor
        #Do a coin flip to determine cis or trans.
        if random.random() < 0.5:
            lastAtom.newCTCenter(otherC, otherN, thisAtom)
        else:
            lastAtom.newCTCenter(otherC, thisAtom, otherN)
    if probablyChiral(lastAtom):
        tempN = []
        for neighbor in lastAtom.neighbors:
            if neighbor != thisAtom:
                tempN.append(neighbor)
        if len(tempN) == 2:
            tempN.append(None)
        #Do a coin flip to determine chirality.
        if random.random() < 0.5:
            lastAtom.newChiralCenter(thisAtom, (tempN[0],
                        tempN[1], tempN[2]))
        else:
            lastAtom.newChiralCenter(thisAtom, (tempN[0],
                        tempN[2], tempN[1]))
        
        
            

def probablyChiral(atom):
    #A rather bootleg heruistic for determining whether an atom
    #may need a chiral center.
    if len(atom.neighbors) <= 2:
        return False
    carbonCount = 0
    elementList = []
    for neighbor in atom.neighbors:
        if neighbor.element == 'C':
            carbonCount += 1
        if neighbor.element not in elementList:
            elementList.append(neighbor.element)
    return carbonCount >= 2 or len(elementList) >= 3

if __name__ == '__main__':
    print smiles(randomStart()[0])

    


