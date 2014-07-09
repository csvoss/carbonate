"""
molecularStructure.py
Contains class Molecule and class Atom.
"""
import copy
##from toMolecule import *

#Testing - replace "H" with "Br" to visualize all hydrogens
hydrogen = "H"

class AlleneError(Exception):
    """
    Raised when we try to make an allene.  Tells higher functions that
    the reaction we just attempted should not be allowed.  (Even if it
    is technically chemically feasible.)
    """
    pass

class Molecule:
    """
    This class represents the structure of a molecule.
        self.atoms :: [Atom].
    """

    def __init__(self, firstAtom):
        """
        firstAtom :: Atom.
        """
        self.atoms = [firstAtom]
    
    def addAtom(self, newAtom, targetAtom, bondOrder=1):
        """
        For adding new atoms to the molecule.

        newAtom :: Atom. The new atom to attach somewhere.
        targetAtom :: Atom. Must be in this molecule. The place newAtom will attach.
        bondOrder :: int. 1 for single, 2 for double, 3 for triple, 4 for quadruple.

        Throws error if targetAtom is not in the molecule.
        Throws error if bondOrder is invalid.
        """
        if bondOrder == 0:
            return
        if bondOrder not in [1,2,3,4]:
            raise StandardError("Invalid bond order: %s" % str(bondOrder))
        if targetAtom not in self.atoms:
            raise StandardError("Error in addAtom: target atom not already in molecule.")
        if newAtom in self.atoms:
            print "WARNING: new atom already in molecule.  Using addBond instead."
            return self.addBond(targetAtom, newAtom, bondOrder)
        self.atoms.append(newAtom)
        self.addBond(newAtom, targetAtom, bondOrder)
        
    def addBond(self, atom1, atom2, bondOrder=1):
        """
        Creates a new bond between two atoms already in the molecule.
        atom1 :: Atom.
        atom2 :: Atom.
        bondOrder :: int. 1, 2, 3, or 4.
        """
        if atom1 not in self.atoms:
            raise StandardError("Error in addBond: atom1 not in molecule. %s" % atom1.element)
        if atom2 not in self.atoms:
            raise StandardError("Error in addBond: atom2 not in molecule. %s" % atom2.element)
        atom1.neighbors[atom2] = bondOrder
        atom2.neighbors[atom1] = bondOrder

    def addMolecule(self, molecule, foreignTarget, selfTarget, bondOrder=1):
        """
        Welds two molecules together by the given atoms.
        molecule :: Molecule.
        foreignTarget :: Atom. The atom on the other molecule.
        selfTarget :: Atom. An atom on this molecule.
        bondOrder :: int. Bond order for the new bond between foreignTarget and selfTarget.
        """
        #Preserves objects in added molecule (no deepcopy)
        for foreignAtom in molecule.atoms:
            if not (foreignAtom in self.atoms):
                self.atoms.append(foreignAtom)
        self.addBond(selfTarget, foreignTarget, bondOrder)
        
    def removeAtom(self, target):
        """
        Remove an atom from this molecule. Destroys the atom.
        target :: Atom.
        """
        for atom in self.atoms:
            if target in atom.neighbors:
                del(atom.neighbors[target])
        self.atoms.remove(target)
        del(target)
        
    def changeBond(self, atom1, atom2, newBondOrder):
        """
        Changes the bond order of a bond.
        atom1 :: Atom.
        atom2 :: Atom.
        newBondOrder :: int. 0, 1, 2, 3, or 4.
        Note that newBondOrder=0 breaks the bond.
        """
        if newBondOrder == 0:
            del(atom1.neighbors[atom2])
            del(atom2.neighbors[atom1])
        else:
            atom1.neighbors[atom2] = newBondOrder
            atom2.neighbors[atom1] = newBondOrder

    def addHydrogens(self):
        """
        Adds a full complement of hydrogens to every atom.
        As a side-effect, also checks for over-valence.
        """
        for atom in self.atoms:
            if atom.element == 'C':
                maxval = 4
            elif atom.element == 'N':
                maxval = 3
            elif atom.element == 'O':
                maxval = 2
            else:
                continue
            val = 0
            for neighbor in atom.neighbors:
                val += atom.neighbors[neighbor]
            if val > maxval:
                print str(atom) + " has too many bonds!"
                raise StandardError
            for i in xrange(maxval-val):
                H = Atom(hydrogen)
                self.addAtom(H, atom, 1)

    def countElement(self, element):
        """
        Counts the occurrences of `element` in this molecule's atoms.
        element :: str. Case matters.
        return :: int.
        """
        out = 0
        for atom in self.atoms:
            if atom.element == element:
                out += 1
        return out

    def removeBond(self, atom1, atom2):
        """
        Remove the bond between atom1 and atom2 without removing either
        from this molecule.
        atom1 :: Atom.
        atom2 :: Atom.
        """
        return self.changeBond(atom1, atom2, 0)

    def withHydrogens(self):
        """
        Return a version of this molecule in which explicit hydrogens have 
        been added to all molecules for which we can infer the number to add.
        """
        BOND_ORDERS = {
            "B": 3,
            "C": 4,
            "N": 3,
            "O": 2,
            "P": 5,
            "S": 2,
            "F": 1,
            "Cl": 1,
            "Br": 1,
            "I": 1,
        }

        out = copy.deepcopy(self)
        for atom in out.atoms:
            elem = atom.element
            bond_order = atom.totalBondOrder()
            charge = atom.charge
            if elem in BOND_ORDERS:
                ideal_bond_order = BOND_ORDERS[elem]
                hydrogens_to_add = ideal_bond_order - bond_order + charge
                for _ in xrange(hydrogens_to_add):
                    new_hydrogen = Atom("H")
                    out.addAtom(new_hydrogen, atom, 1)
        return out



class Atom:
    """
    This class represents the structure of a single atom. Not to be used alone;
    only really useful as a part of a Molecule.
        self.element :: str.
        self.charge :: int.
        self.neighbors :: {Atom, int}.  This atom's neighbors; the int represents bond order.

    Refer to moleculeToSmiles for these attributes' purpose; **ignore** them otherwise.
        self.flag :: int.
        self.rflag :: [(int, Atom)].
        self.nRead :: int.
        self.parentAtom :: 0 or Atom.
        self.nonHNeighbors = [Atom].
    """

    def __str__(self):

        ## Atoms are represented by the standard abbreviation of the chemical elements, in square brackets, such as [Au] for gold. Brackets can be omitted for the "organic subset" of B, C, N, O, P, S, F, Cl, Br, and I. All other elements must be enclosed in brackets. If the brackets are omitted, the proper number of implicit hydrogen atoms is assumed; for instance the SMILES for water is simply O.

        ## An atom holding one or more electrical charges is enclosed in brackets, followed by the symbol H if it is bonded to one or more atoms of hydrogen, followed by the number of hydrogen atoms (as usual one is omitted example: NH4 for ammonium), then by the sign '+' for a positive charge or by '-' for a negative charge. The number of charges is specified after the sign (except if there is one only); however, it is also possible write the sign as many times as the ion has charges: instead of "Ti+4", one can also write "Ti++++" (Titanium IV, Ti4+). Thus, the hydroxide anion is represented by [OH-], the oxonium cation is [OH3+] and the cobalt III cation (Co3+) is either [Co+3] or [Co+++].

        ## ALL MOLECULES SHOULD HAVE HYDROGENS EXPLICITLY ATTACHED AS ATOMS

        ## We're ignoring isotopes.

        if startAtom.charge < 0:
            outp = "["+outp + str(startAtom.charge)+"]"
        if startAtom.charge > 0:
            outp = "["+outp + "+"+str(startAtom.charge)+"]"

        return self.element
    
    def __init__(self, element, charge=0):
        """
        Initialize. Charge, by default, is 0.
        element :: str.  Should be the periodic table abbreviation.
        charge :: int. Optional.
        """
        self.element = element
        self.charge = 0
        self.neighbors = dict()

        #Temporary values which should only be meaningful within smiles() and subsmiles().
        #for traversing
        self.flag = 0
        #for ring-finding
        self.rflag = [] #empty if not part of a ring bond, nonempty if otherwise
        self.nRead = 0 #neighbors already read
        self.parentAtom = 0 #atom right before this one
        self.nonHNeighbors = []

    def newChiralCenter(self, reference, clockwiseList):
        """
        Set up this atom as a chiral center.
        reference :: Atom.
        clockwiseList :: a list of 3 Atoms.
        """
        self.chiralA = reference
        self.chiralB, self.chiralC, self.chiralD = clockwiseList

    def chiralCWlist(self, reference):
        """
        Returns a list of the other 3 Atoms bonded to this Atom,
        in **clockwise** order when looking down `reference`.
        Note: Left-hand rule instead of right-hand rule.

        reference :: Atom.
        return :: a list of 3 Atoms.
        """
        if reference == self.chiralA:
            return [self.chiralB, self.chiralC, self.chiralD]
        elif reference == self.chiralB:
            return [self.chiralA, self.chiralD, self.chiralC]
        elif reference == self.chiralC:
            return [self.chiralA, self.chiralB, self.chiralD]
        elif reference == self.chiralD:
            return [self.chiralA, self.chiralC, self.chiralB]
        else:
            msg = "Error in chiralCWlist: no such reference atom: %s\n" % reference.element
            msg += ", ".join[self.chiralA.element, self.chiralB.element,
                             self.chiralC.element, self.chiralD.element]
            raise StandardError(msg)
    
    def chiralRingList(self, inport, outport):
        """
        Returns which substituent is up, followed by which one is down,
        in a ring context. Assume that the ring veers left (ccw).

        inport :: Atom.
        outport :: Atom.
        return :: a tuple of 2 Atoms.

        **DEPRECATED**. Please stick to using chiralCWlist.
        """
        if inport == self.chiralA:
            if outport == self.chiralB:
                return                  (self.chiralC, self.chiralD)
            if outport == self.chiralC:
                return                  (self.chiralD, self.chiralB)
            if outport == self.chiralD:
                return                  (self.chiralB, self.chiralC)
        elif inport == self.chiralB:
            if outport == self.chiralA:
                return                  (self.chiralD, self.chiralC)
            if outport == self.chiralC:
                return                  (self.chiralA, self.chiralD)
            if outport == self.chiralD:
                return                  (self.chiralC, self.chiralA)
        elif inport == self.chiralC:
            if outport == self.chiralA:
                return                  (self.chiralB, self.chiralD)
            if outport == self.chiralB:
                return                  (self.chiralD, self.chiralA)
            if outport == self.chiralD:
                return                  (self.chiralA, self.chiralB)
        elif inport == self.chiralD:
            if outport == self.chiralA:
                return                  (self.chiralC, self.chiralB)
            if outport == self.chiralB:
                return                  (self.chiralA, self.chiralC)
            if outport == self.chiralC:
                return                  (self.chiralB, self.chiralA)
                
        msg = "Error in chiralCWlist: no such inport and outport: %s, %s\n" %\
              (inport.element, outport.element)
        msg += ", ".join[self.chiralA.element, self.chiralB.element,
                         self.chiralC.element, self.chiralD.element]
        raise StandardError(msg)
        
    def newCTCenter(self, otherC, a, b):
        """
        CTCenters (cis-trans centers) must come in pairs.  Both of the
        carbons across the double bond must have a CTCenter -- so you must
        execute this method once for each. Assuming that you're using the same
        plane of reference for each CT-center, this makes Atom a to be
        directly clockwise from otherC, and Atom b directly counterclockwise.

        otherC :: Atom.
        a :: Atom.
        b :: Atom.

        Raises AlleneError if newCTCenter has already been applied to this molecule.
        """
        if hasattr(self, 'CTa'):
            ## That would make this carbon the center of an allene!
            ## Keeping track of stereochem for allenes is Hard.
            ## TODO: Make it so that we can in fact handle allenes?
            raise AlleneError
        self.CTotherC = otherC
        self.CTa = a
        self.CTb = b

        if not (self.neighbors[otherC] == 2 and otherC.neighbors[self] == 2):
            raise StandardError("Error in newCTCenter: cis-trans center without double bond")
    
    def eliminateChiral(self):
        """
        Destroys the chirality information. Don't worry, the atoms are still there.
        """
        if hasattr(self, 'chiralA'):
            del(self.chiralA)
            del(self.chiralB)
            del(self.chiralC)
            del(self.chiralD)

    def eliminateCT(self):
        """
        Destroys the cis-trans information. Don't worry, the atoms are still there.
        """
        del(self.CTotherC)
        del(self.CTa)
        del(self.CTb)

    def totalBondOrder(self):
        """
        Counts neighbors and returns the total bond order of this atom,
        not including implicit hydrogens.

        return :: int.
        """
        out = 0
        for neighbor in self.neighbors:
            out += self.neighbors[neighbor]
        return out

    def findAlkeneBond(self):
        """
        Returns the carbon atom to which this one is double-bonded, or None if no such.

        return :: Atom or None.
        """
        for neighbor in self.neighbors:
            if self.neighbors[neighbor] == 2 and neighbor.element == 'C':
                return neighbor
        return None
        



## Example usage
"""
c40 = Atom("C")
c41 = Atom("C")
mol4 = Molecule(c40)
mol4.addAtom(c41, c40, 2)
c42 = Atom("C")
c43 = Atom("C")
mol4.addAtom(c42, c40, 1)
mol4.addAtom(c43, c41, 1)
c40.newCTCenter(c41, c42, None)
c41.newCTCenter(c40, c43, None)
"""