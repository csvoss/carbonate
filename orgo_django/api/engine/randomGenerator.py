from helperFunctions import *
import random
from toMolecule import moleculify
from toSmiles import smilesify

def random_molecule(end_prob=0.3, max_length=15, alkyne_prob=0.2,
                    alkeneProb=0.4, br_prob=0.1, cl_prob=0.1, oh_prob=0.1,
                    branch_prob=0.05, force_terminal_alkyne=False):
    mol, _, _ = random_start(end_prob, max_length, alkyne_prob,
                             alkeneProb, br_prob, cl_prob, oh_prob,
                             branch_prob, force_terminal_alkyne)
    print "Made the molecule"
    mol.addHydrogens()
    print "Added hydrogens"
    return mol

def random_start(end_prob=0.3, max_length=15, alkyne_prob=0.2,
                 alkeneProb=0.4, br_prob=0.1, cl_prob=0.1, oh_prob=0.1,
                 branch_prob=0.05, force_terminal_alkyne=False):
    #Want more alkynes?  You need to make fewer substituents.
    print "random_start"
    for prob in (br_prob, cl_prob, oh_prob, branch_prob):
        prob *= 1 - 2*alkyne_prob
    last_atom = Atom("C")
    front_atom = last_atom
    mol = Molecule(last_atom)
    if force_terminal_alkyne:
        last_atom = Atom("C")
        mol.addAtom(last_atom, front_atom, 3)
    while (random.random() < 1.0-end_prob or len(mol.atoms) < 3) and len(mol.atoms) < max_length:
        switcher = random.random()
        if switcher < 0.8:
            newMol, atom, nextAtom = randC(br_prob, cl_prob, oh_prob, branch_prob)
        elif switcher < 0.9:
            newMol, atom, nextAtom = random_ring(5, br_prob, cl_prob, oh_prob)
        else:
            newMol, atom, nextAtom = random_ring(6, br_prob, cl_prob, oh_prob)

        switcher = random.random()
        if switcher < alkyne_prob and len(atom.neighbors) == 0\
           and last_atom.totalBondOrder() == 1:
            #Make a triple bond.
            mol.addMolecule(newMol, atom, last_atom, 3)
        elif switcher < alkeneProb+alkyne_prob and len(atom.neighbors) <= 1\
             and last_atom.totalBondOrder() <= 2\
             and last_atom.findAlkeneBond() == None\
             and last_atom.totalBondOrder() > 0:
            #Make a double bond.  Watch out for cis/trans.
            #Each double bond must have a cis/trans specification.
            mol.addMolecule(newMol, atom, last_atom, 2)
            others = []
            for neighbor in last_atom.neighbors:
                if neighbor != atom:
                    others.append(neighbor)
            if len(others) == 1:
                last_atom.newCTCenter(atom, others[0], Atom("H"))
            else:
                last_atom.newCTCenter(atom, others[0], others[1])
            #Flag atom so that we add stereo to it, later.
            atom.makeCTFlag = True
        else:
            #Make single bond
            mol.addMolecule(newMol, atom, last_atom, 1)
        fix_stereo(mol, atom, last_atom)
        last_atom = atom
        atom = nextAtom
    #Not quite done yet.  We need to fix stereochemistry on the last atom.
    fix_stereo(mol, Atom("H"), last_atom) #May be slight problem with chirality?
    return mol, front_atom, None

def randC(cl_prob, br_prob, oh_prob, branch_prob):
    '''
    Makes a random carbon atom with some substituents
    '''
    c = Atom("C")
    mol = Molecule(c)
    #Add up to 2 substituents
    for i in xrange(2):
        switcher = random.random()
        if switcher < cl_prob:
            newS = Atom("Cl")
            mol.addAtom(newS, c, 1)
        elif switcher < cl_prob + br_prob:
            newS = Atom("Br")
            mol.addAtom(newS, c, 1)
        elif switcher < cl_prob + br_prob + oh_prob:
            newS = Atom("O")
            mol.addAtom(newS, c, 1)
        elif switcher < cl_prob + br_prob + oh_prob + branch_prob:
            print "Branching"
            newS, front_atom, _ = random_start()
            mol.addMolecule(newS, front_atom, c, 1)
    return mol, c, c

def random_ring(number_of_carbons, br_prob, cl_prob, oh_prob):
    '''
    Makes a random ring
    '''
    init_atom = Atom("C")
    mol = Molecule(init_atom)
    old_atom = init_atom
    for _ in xrange(number_of_carbons - 1):
        new_atom = Atom('C')
        mol.addAtom(new_atom, old_atom, 1)
        old_atom = new_atom
    #Now, close the loop!
    mol.addBond(init_atom, old_atom, 1)
    out_atom = mol.atoms[random.randint(1, len(mol.atoms)-1)]
    for atom in mol.atoms:
        #Don't interfere with the larger structure.
        if atom == init_atom or atom == out_atom:
            continue
        #Add some random substituents.
        switcher = random.random()
        if switcher < br_prob:
            addAtom = Atom("Br")
        elif switcher < br_prob + cl_prob:
            addAtom = Atom("Cl")
        elif switcher < br_prob + cl_prob + oh_prob:
            addAtom = Atom("O")
        else:
            continue
        if len(atom.neighbors) < 2:
            continue
        #Do a coin flip to determine chirality.
        if random.random() < 0.5:
            atom.newChiralCenter(addAtom, (atom.neighbors.keys()[0], atom.neighbors.keys()[1], Atom("H")))
        else:
            atom.newChiralCenter(addAtom, (atom.neighbors.keys()[1], atom.neighbors.keys()[0], Atom("H")))
        mol.addAtom(addAtom, atom, 1)
    
    return mol, init_atom, out_atom

def fix_stereo(mol, atom, last_atom):
    '''
    Look to see if the *last* piece we added requires stereochem
    '''
    otherC = last_atom.findAlkeneBond()
    if hasattr(last_atom, "makeCTFlag"):
        del last_atom.makeCTFlag
        #Find the last neighbor of last_atom (other than atom
        #and otherC).  None means Hydrogen.
        otherN = Atom("H")
        for neighbor in last_atom.neighbors:
            if neighbor != otherC and neighbor != atom:
                otherN = neighbor
        #Do a coin flip to determine cis or trans.
        if random.random() < 0.5:
            last_atom.newCTCenter(otherC, otherN, atom)
        else:
            last_atom.newCTCenter(otherC, atom, otherN)
    if probablyChiral(last_atom):
        tempN = []
        for neighbor in last_atom.neighbors:
            if neighbor != atom:
                tempN.append(neighbor)
        if len(tempN) == 2:
            new_h = Atom("H")
            mol.addAtom(new_h, last_atom, 1)
            tempN.append(new_h)
        #Do a coin flip to determine chirality.
        if random.random() < 0.5:
            last_atom.newChiralCenter(atom, (tempN[0], tempN[1], tempN[2]))
        else:
            last_atom.newChiralCenter(atom, (tempN[0], tempN[2], tempN[1]))
        
        
            

def probablyChiral(atom):
    '''
    A rather bootleg heuristic for determining whether an atom is chiral.
    '''
    #may need a chiral center.
    if len(atom.neighbors) <= 2:
        return False
    carbon_count = 0
    element_list = []
    for neighbor in atom.neighbors:
        if neighbor.element == 'C':
            carbon_count += 1
        if neighbor.element not in element_list:
            element_list.append(neighbor.element)
    return carbon_count >= 2 or len(element_list) >= 3



