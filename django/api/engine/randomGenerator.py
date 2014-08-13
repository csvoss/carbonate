import random

from molecularStructure import Atom, Molecule

ALKANEPROB = 1.0
terminate_count = 0

def update_terminate_count():
    global terminate_count
    updater = random.random()
    diff = 1 - terminate_count
    terminate_count = 1 - updater * diff

def reset_terminate_count():
    global terminate_count
    terminate_count = 0

def get_terminate_count():
    global terminate_count
    return terminate_count

def random_molecule(force_terminal_alkyne=False, seed=None):

    # Only use seed for debugging
    # if seed is not None:
    #     random.seed(seed)

    global carboncount
    carboncount = 0

    start_atom = Atom("C")
    start_atom.clss = carboncount
    carboncount += 1

    molecule = Molecule(start_atom)

    if force_terminal_alkyne:
        alkyne_atom = Atom("C")
        alkyne_atom.clss = carboncount
        carboncount += 1
        molecule.addAtom(alkyne_atom, start_atom, 3)
        start_atom = alkyne_atom

    global brprob, clprob, ohprob, hprob, nprob
    global etherprob, amineprob, alkeneprob, alkyneprob
    alkeneprob = 0.5 * random.random()
    alkyneprob = 0.2 * random.random()
    brprob = 0.1 * random.random()
    clprob = 0.1 * random.random()
    ohprob = 0.1 * random.random()
    nprob = 0.1 * random.random()
    hprob = 0.5 * random.random()
    etherprob = 0.1 * random.random()
    amineprob = 0.05 * random.random()

    tree, end_atom = random_tree()
    molecule.addMolecule(tree, end_atom, start_atom, 1)

    molecule.addHydrogens()

    add_random_links(molecule)
    add_random_chirality(molecule)
    add_random_cistrans(molecule)

    reset_terminate_count()

    return molecule

def random_tree(ether=True, amine=True):
    """
    Precondition: The atom created herein will only have one bond to it
    to its parent atom.
    return :: Molecule, Atom.
    """
    global etherprob, amineprob, carboncount
    rand = random.random()
    if rand < etherprob and ether:
        atom = Atom("O")
        molecule = Molecule(atom)
        new_molecule, new_atom = random_tree(ether=False, amine=False)
        molecule.addMolecule(new_molecule, new_atom, atom, 1)
        return molecule, atom
    elif rand < etherprob+amineprob and amine:
        atom = Atom("N")
        molecule = Molecule(atom)
        new_molecule, new_atom = random_tree(ether=False, amine=False)
        molecule.addMolecule(new_molecule, new_atom, atom, 1)
        new_molecule, new_atom = random_tree(ether=False, amine=False)
        molecule.addMolecule(new_molecule, new_atom, atom, 1)
        return molecule, atom

    atom = Atom("C")
    atom.clss = carboncount
    carboncount += 1
    molecule = Molecule(atom)

    global alkeneprob, alkyneprob
    rand = random.random()
    totalprob = float(alkeneprob + alkyneprob + ALKANEPROB)
    if rand < alkeneprob / totalprob:
        bond_orders = [2, 1]
    elif rand < (alkyneprob + alkeneprob) / totalprob:
        bond_orders = [3]
    else:
        bond_orders = [1, 1, 1]

    for bond_order in bond_orders:
        new_molecule, new_atom = random_constituent(bond_order)
        molecule.addMolecule(new_molecule, new_atom, atom, bond_order)

    return molecule, atom

def random_constituent(n):
    """
    Return a molecular branch to attach to a parent atom.
    Compatible with the bond to parent atom being of order n.
    return :: Molecule, Atom
    """
    global carboncount
    update_terminate_count()
    rand = random.random()
    if rand < get_terminate_count():
        if n == 1:
            atom = random_terminal()
        elif n == 2 or n == 3:
            atom = Atom("C")
            atom.clss = carboncount
            carboncount += 1
        else:
            raise StandardError("Invalid number of constituents to add")
        return Molecule(atom), atom
    else:
        return random_tree_safe()

def random_terminal():
    """
    Return a random single atom (a "terminal"). No branching.
    return :: Atom
    """
    global ohprob, brprob, clprob, hprob, nprob
    rand = random.random()
    totalprob = float(ohprob + brprob + clprob + hprob + nprob)
    if rand < ohprob / totalprob:
        return Atom("O")
    elif rand < (ohprob + brprob) / totalprob:
        return Atom("Br")
    elif rand < (ohprob + brprob + clprob) / totalprob:
        return Atom("Cl")
    elif rand < (ohprob + brprob + clprob + nprob) / totalprob:
        return Atom("N")
    else:
        return Atom("H")

def random_tree_safe():
    """
    Return a molecular branch in which the first atom has only a single
    bond to another carbon -- safe for attachment to alkynes (and alkenes, too).
    TODO: Make a variant of this which is more suited for alkenes and will
          allow a greater diversity of substituents attached to alkenes.
    return :: Molecule, Atom.
    """
    global carboncount
    atom = Atom("C")
    atom.clss = carboncount
    carboncount += 1
    molecule = Molecule(atom)
    atom2 = Atom("C")
    atom2.clss = carboncount
    carboncount += 1
    molecule.addAtom(atom2, atom, 1)
    new_molecule, new_atom = random_tree()
    molecule.addMolecule(new_molecule, new_atom, atom2, 1)
    return molecule, atom

def add_random_links(molecule):
    """
    Add a ring if one is possible to add.
    return :: None.
    """
    ring_size = random.choice([5, 6])

    carbons = [i for i in molecule.atoms if i.element is 'C']
    carbons = shuffled(carbons)

    def get_n_away(parent, start_carbon, n):
        noncarbons = [i for i in start_carbon.neighbors if i.element != 'C']
        if n == 0:
            if len(noncarbons) > 0 and start_carbon.element == 'C':
                return start_carbon, random.choice(noncarbons)
            else: ## self does not have a noncarbon substituent, or is not carbon
                return None, None
        for next_atom in shuffled(start_carbon.neighbors):
            if next_atom in carbons and next_atom is not parent:
                if start_carbon.neighbors[next_atom] < 3:
                    return get_n_away(start_carbon, next_atom, n-1)
        return None, None

    ## Add a ring at random if it is discovered to be possible
    if len(carbons) > 4:
        for start_carbon in carbons:
            noncarbons = [i for i in start_carbon.neighbors if i.element != 'C']
            if len(noncarbons) <= 0:
                continue
            start_noncarbon = random.choice(noncarbons)
            for next_atom in shuffled(start_carbon.neighbors):
                if next_atom in carbons:
                    if start_carbon.neighbors[next_atom] < 3:
                        result_carbon, result_noncarbon = get_n_away(start_carbon, next_atom, ring_size-2)
                        if result_carbon == None:
                            return
                        else:
                            ## delete a random noncarbon substituent from each atom
                            molecule.removeAtom(start_noncarbon)
                            molecule.removeAtom(result_noncarbon)
                            ## make a bond between the two atoms
                            molecule.addBond(start_carbon, result_carbon, 1)
                            return

def add_random_chirality(molecule):
    """
    return :: None.
    """
    ## MUST BE DONE **AFTER** ADDING LINKS
    ## TODO
    pass

def add_random_cistrans(molecule):
    """
    Add arbitrary cis/trans centers where relevant.
    return :: None.
    """
    for atom in molecule.atoms:
        for other_atom in atom.neighbors:
            if atom.neighbors[other_atom] == 2:
                atoms = [a for a in atom.neighbors if a is not other_atom]
                if len(atoms) == 2:
                    atom.newCTCenter(other_atom, atoms[0], atoms[1])
                else:
                    raise StandardError("Too many constituents on = carbon")
                atoms = [a for a in other_atom.neighbors if a is not atom]
                if len(atoms) == 2:
                    other_atom.newCTCenter(atom, atoms[0], atoms[1])
                else:
                    raise StandardError("Too many constituents on = carbon")

def shuffled(l):
    """
    return :: list.
    """
    new_l = [i for i in l]
    if len(new_l) > 1:
        random.shuffle(new_l)
    return new_l