# import random

from toMolecule import moleculify
from toSmiles import smilesify
from molecularStructure import Atom, Molecule

random_random = lambda: 0.8 ## TEST
random_randint = lambda foo, bar: 1 ## TEST

carbon_counter = 0

def random_molecule(end_prob=0.3, max_length=15, alkyne_prob=0.2,
                    alkeneProb=0.4, br_prob=0.1, cl_prob=0.1, oh_prob=0.1,
                    branch_prob=0.05, force_terminal_alkyne=False):
    mol, _, _ = random_start(end_prob, max_length, alkyne_prob,
                             alkeneProb, br_prob, cl_prob, oh_prob,
                             branch_prob, force_terminal_alkyne)
    mol.addHydrogens()
    return mol

def random_start(end_prob=0.3, max_length=15, alkyne_prob=0.2,
                 alkeneProb=0.4, br_prob=0.1, cl_prob=0.1, oh_prob=0.1,
                 branch_prob=0.05, force_terminal_alkyne=False):
    global carbon_counter
    #Want more alkynes?  You need to make fewer substituents.
    # for prob in (br_prob, cl_prob, oh_prob, branch_prob):
    #     prob *= 1 - 2*alkyne_prob
    last_atom = Atom("C")
    last_atom.clss = carbon_counter
    carbon_counter += 1
    front_atom = last_atom
    mol = Molecule(last_atom)
    if force_terminal_alkyne:
        last_atom = Atom("C")
        last_atom.clss = carbon_counter
        carbon_counter += 1
        mol.addAtom(last_atom, front_atom, 3)
    while (random_random() < 1.0-end_prob or len(mol.atoms) < 3) and len(mol.atoms) < max_length:
        switcher = random_random()
        if switcher < 0.8:
            newMol, atom, nextAtom = rand_carbon(br_prob, cl_prob, oh_prob, branch_prob)
        elif switcher < 0.9:
            newMol, atom, nextAtom = random_ring(5, br_prob, cl_prob, oh_prob)
        else:
            newMol, atom, nextAtom = random_ring(6, br_prob, cl_prob, oh_prob)

        switcher = random_random()
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

def rand_carbon(cl_prob, br_prob, oh_prob, branch_prob):
    '''
    Makes a random carbon atom with some substituents
    '''
    global carbon_counter
    carbon = Atom("C")
    carbon.clss = carbon_counter
    carbon_counter += 1
    mol = Molecule(carbon)
    #Add up to 3 substituents
    for _ in xrange(3):
        switcher = random_random()
        if switcher < cl_prob:
            next_branch = Atom("Cl")
            mol.addAtom(next_branch, carbon, 1)
        elif switcher < cl_prob + br_prob:
            next_branch = Atom("Br")
            mol.addAtom(next_branch, carbon, 1)
        elif switcher < cl_prob + br_prob + oh_prob:
            next_branch = Atom("O")
            mol.addAtom(next_branch, carbon, 1)
        elif switcher < cl_prob + br_prob + oh_prob + branch_prob:
            next_branch, first_atom, _ = random_start()
            mol.addMolecule(next_branch, first_atom, carbon, 1)
            mol.addAtom(Atom("Cf"), first_atom, 1) ## TEST
            mol.addAtom(Atom("Es"), carbon, 1) ## TEST
        else:
            # mol.addAtom(Atom("H"), carbon, 1)
            pass
    return mol, carbon, carbon

def random_ring(number_of_carbons, br_prob, cl_prob, oh_prob):
    '''
    Makes a random ring
    '''
    global carbon_counter
    init_atom = Atom("C")
    init_atom.clss = carbon_counter
    carbon_counter += 1
    mol = Molecule(init_atom)
    old_atom = init_atom
    for _ in xrange(number_of_carbons - 1):
        new_atom = Atom("C")
        new_atom.clss = carbon_counter
        carbon_counter += 1
        mol.addAtom(new_atom, old_atom, 1)
        old_atom = new_atom
    #Now, close the loop!
    mol.addBond(init_atom, old_atom, 1)
    rand_index = random_randint(1, len(mol.atoms)-1)
    out_atom = mol.atoms[rand_index]
    for atom in mol.atoms:
        #Don't interfere with the larger structure.
        if atom is init_atom or atom is out_atom:
            continue
        #Add some random substituents.
        switcher = random_random()
        if switcher < br_prob:
            new_atom = Atom("Br")
        elif switcher < br_prob + cl_prob:
            new_atom = Atom("Cl")
        elif switcher < br_prob + cl_prob + oh_prob:
            new_atom = Atom("O")
        else:
            new_atom = Atom("H") ## before: was a continue
        if len(atom.neighbors) < 3:
            continue
        #Do a coin flip to determine chirality.
        if random_random() < 0.5:
            atom.newChiralCenter(new_atom, (atom.neighbors.keys()[0], atom.neighbors.keys()[1], Atom("H")))
        else:
            atom.newChiralCenter(new_atom, (atom.neighbors.keys()[1], atom.neighbors.keys()[0], Atom("H")))

        mol.addAtom(new_atom, atom, 1)
    
    return mol, init_atom, out_atom

def fix_stereo(mol, atom, last_atom):
    '''
    Look to see if the *last* piece we added requires stereochem
    '''
    other_carbon = last_atom.findAlkeneBond()
    if hasattr(last_atom, "makeCTFlag"):
        del last_atom.makeCTFlag
        #Find the last neighbor of last_atom (other than atom
        #and other_carbon).  None means Hydrogen.
        other_neighbor = Atom("H")
        for neighbor in last_atom.neighbors:
            if neighbor != other_carbon and neighbor != atom:
                other_neighbor = neighbor
        #Do a coin flip to determine cis or trans.
        if random_random() < 0.5:
            last_atom.newCTCenter(other_carbon, other_neighbor, atom)
        else:
            last_atom.newCTCenter(other_carbon, atom, other_neighbor)
    
    elif probably_chiral(last_atom):
        temp_neighbors = []
        for neighbor in last_atom.neighbors:
            if neighbor is not atom:
                temp_neighbors.append(neighbor)
        if len(temp_neighbors) == 2:
            new_h = Atom("H")
            mol.addAtom(new_h, last_atom, 1)
            mol.addAtom(Atom("Cu"), new_h, 1) ## TEST
            temp_neighbors.append(new_h)
        #Do a coin flip to determine chirality.
        if random_random() < 0.5:
            last_atom.newChiralCenter(atom, (temp_neighbors[0], temp_neighbors[1], temp_neighbors[2]))
        else:
            last_atom.newChiralCenter(atom, (temp_neighbors[0], temp_neighbors[2], temp_neighbors[1]))


def probably_chiral(atom):
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
