#Testing - replace "H" with "Br" to visualize all hydrogens
hydrogen = "H"

debugSmiles = False

class AlleneError(Exception):
    #Raised when we try to make an allene.  Tells higher functions that
    #the reaction we just attempted should not be allowed.  (Even if it
    #is technically chemically feasible.)
    pass

class Molecule:
    
    def __init__(self, firstAtom):
        self.atoms = [firstAtom]

    
    def addAtom(self, newAtom, targetAtom, bondOrder):
        if targetAtom not in self.atoms:
            print "Error in addAtom: target atom not already in molecule."
            raise StandardError
        if newAtom in self.atoms:
            print "Error in addAtom: new atom already in molecule.  Use addBond instead."
            raise StandardError
        self.atoms.append(newAtom)
        self.addBond(newAtom, targetAtom, bondOrder)
        
    def addBond(self, atom1, atom2, bondOrder):
        atom1.neighbors[atom2] = bondOrder
        atom2.neighbors[atom1] = bondOrder

    def addMolecule(self, molecule, foreignTarget, selfTarget, bo):
        #Preserves objects in added molecule (no deepcopy)
        for foreignAtom in molecule.atoms:
            if not (foreignAtom in self.atoms):
                self.atoms.append(foreignAtom)
        self.addBond(selfTarget, foreignTarget, bo)
        
    def removeAtom(self, target):
        for atom in self.atoms:
            if target in atom.neighbors:
                del(atom.neighbors[target])
        self.atoms.remove(target)
        del(target)
        
    def changeBond(self, atom1, atom2, newBondOrder):
        #newBondOrder=0 breaks the bond.
        if newBondOrder == 0:
            del(atom1.neighbors[atom2])
            del(atom2.neighbors[atom1])
        else:
            atom1.neighbors[atom2] = newBondOrder
            atom2.neighbors[atom1] = newBondOrder

    def addHydrogens(self):
        #Adds a full complement of hydrogens to every atom.
        #As a side-effect, also checks for over-valence.
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
        out = 0
        for atom in self.atoms:
            if atom.element == element:
                out += 1
        return out


class Atom:

##    def __str__(self):
##        return self.element
    
    def __init__(self, element):
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
        #Set up this atom as a chiral center.
        #reference is an Atom; clockwiseList is a list of 3 Atoms
        self.chiralA = reference
        self.chiralB, self.chiralC, self.chiralD = clockwiseList

    def chiralCWlist(self, reference):
        #Returns a list of the other 3 Atoms bonded to this Atom,
        #in clockwise order when looking down reference.
        if reference == self.chiralA:
            return [self.chiralB, self.chiralC, self.chiralD]
        elif reference == self.chiralB:
            return [self.chiralA, self.chiralD, self.chiralC]
        elif reference == self.chiralC:
            return [self.chiralA, self.chiralB, self.chiralD]
        elif reference == self.chiralD:
            return [self.chiralA, self.chiralC, self.chiralB]
        else:
            print "chiralCWlist: no such reference."
            print reference.element
            print reference
            print self.chiralA, self.chiralB, self.chiralC, self.chiralD
            raise StandardError
    
    def chiralRingList(self, inport, outport):
        #Returns which substituent is up, followed by which one is down,
        #in a ring context.
        if inport == self.chiralA:
            if outport == self.chiralB:
                return                  (self.chiralC, self.chiralD)
            if outport == self.chiralC:
                return                  (self.chiralB, self.chiralD)
            if outport == self.chiralD:
                return                  (self.chiralB, self.chiralC)
        elif inport == self.chiralB:
            if outport == self.chiralA:
                return                  (self.chiralD, self.chiralC)
            if outport == self.chiralC:
                return                  (self.chiralA, self.chiralD)
            if outport == self.chiralD:
                return                  (self.chiralA, self.chiralC)
        elif inport == self.chiralC:
            if outport == self.chiralA:
                return                  (self.chiralD, self.chiralB)
            if outport == self.chiralB:
                return                  (self.chiralD, self.chiralA)
            if outport == self.chiralD:
                return                  (self.chiralB, self.chiralA)
        elif inport == self.chiralD:
            if outport == self.chiralA:
                return                  (self.chiralC, self.chiralB)
            if outport == self.chiralB:
                return                  (self.chiralC, self.chiralA)
            if outport == self.chiralC:
                return                  (self.chiralA, self.chiralB)
                
        raise StandardError
        
    def newCTCenter(self, otherC, a, b):
        #CTCenters (cis-trans centers) must come in pairs.  Both of the
        #carbons across the double bond must have a CTCenter.  Atom a is
        #directly clockwise from otherC.  Atom b is directly counterclockwise.
        if hasattr(self, 'CTa'):
            print "-----allene error-----"
            raise AlleneError
        self.CTotherC = otherC
        self.CTa = a
        self.CTb = b
    
    def eliminateChiral(self):
        if hasattr(self, 'chiralA'):
            del(self.chiralA)
            del(self.chiralB)
            del(self.chiralC)
            del(self.chiralD)

    def eliminateCT(self):
        del(self.CTotherC)
        del(self.CTa)
        del(self.CTb)

    def totalBondOrder(self):
        #Returns the total bond order, not including hydrogens
        out = 0
        for neighbor in self.neighbors:
            out += self.neighbors[neighbor]
        return out

    def findAlkeneBond(self):
        #Returns the carbon atom to which this one is double-bonded.
        for neighbor in self.neighbors:
            if self.neighbors[neighbor] == 2 and neighbor.element == 'C':
                return neighbor
        return None
        



def smiles(molecule):

    if isinstance(molecule, list):
        return [smiles(molec) for molec in molecule]

    
    if len(molecule.atoms)==0: return ""
    
    ringsfound = 0
    curAtom = molecule.atoms[0]
    homeAtom = molecule.atoms[0]

    #Create the dictionary nonHNeighbors for each atom.
    for atom in molecule.atoms:
        atom.nonHNeighbors = dict((a, b) for (a, b)in atom.neighbors.items() if a.element.lower() != "h")
    

    #Traverse the molecule once, to hunt down and flag rings.
    #Each iteration: (...while we aren't back to the home atom, or if we are,
                    #while the home atom still has neighbors to read)
    while ((curAtom != homeAtom) or (homeAtom.nRead < len(homeAtom.nonHNeighbors))):
    
        #flag current atom as "read" (flag = 1)
        curAtom.flag = 1
        if debugSmiles:
            print "Flagged "+str(curAtom)+" as read."
        
        #if there are neighbors left to read from this atom:
        if (curAtom.nRead < len(curAtom.nonHNeighbors)):
            
            #if the next atom is the parent atom:
            if list(curAtom.nonHNeighbors)[curAtom.nRead] == curAtom.parentAtom:
                #don't do anything but incrementing nRead
                curAtom.nRead += 1
                if debugSmiles:
                    print "Nope, not progressing to parent."
                
            #else,
            else:
                if list(curAtom.nonHNeighbors)[curAtom.nRead].nRead == 0:
                #if the next atom has not been traversed already:
                    curAtom.nRead += 1
                    list(curAtom.nonHNeighbors)[curAtom.nRead - 1].parentAtom = curAtom
                    curAtom = list(curAtom.nonHNeighbors)[curAtom.nRead - 1]
                    #increment current atom's nRead counter
                    #make the next atom the current atom:
                    #make the old atom the next atom's parent
                    if debugSmiles:
                        print "Progressing..."

                else:
                #if the next atom has been traversed already:
                    #it's a ring!
                    ringsfound += 1
                    curAtom.rflag += [(ringsfound, list(curAtom.nonHNeighbors)[curAtom.nRead])]
                    list(curAtom.nonHNeighbors)[curAtom.nRead].rflag += [(ringsfound, curAtom)]
                    curAtom.nRead += 1
                    #increment ringsfound
                    #set rflag on both atoms to ringsfound
                    #increment current atom's nRead counter
                    #make sure the atoms know who each other is
                    if debugSmiles:
                        print "Ring "+str(ringsfound)+" found! "+str(curAtom.rflag)+", "+str(list(curAtom.nonHNeighbors)[curAtom.nRead - 1])
                
        #if not:
            #go backwards to parent atom:
            #set curAtom to its parent atom
        else:
            if debugSmiles:
                print "Regressing to "+str(curAtom.parentAtom)+" from "+str(curAtom)
            curAtom = curAtom.parentAtom
            
            
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

    return outp

bondSymbols = ['0', '-', '=', '#', '4', '5', '6', '7', '8', '9']

#Precondition: molecule has been flagged for ring positioning (some rflag values on atoms might != 0). This is done by smiles().
#Creates and returns a SMILES string for unflagged (!atom.flag==2) atoms within a molecule, starting with the given atom.
def subsmiles(molecule, startAtom, parentAtom):

    if debugSmiles:
        if parentAtom != 0:
            print "\t\t\t [ Entering subsmiles with "+str((startAtom.element, startAtom, startAtom.rflag))+" from "+str((parentAtom.element, parentAtom, parentAtom.rflag))
        else:
            print "\t\t\t [ Entering subsmiles with "+str((startAtom.element, startAtom, startAtom.rflag))+"..."
    
    #Flag the current atom.
    startAtom.flag = 2

    outp = startAtom.element

    if startAtom.charge < 0:
        outp = "["+outp + str(startAtom.charge)+"]"
    if startAtom.charge > 0:
        outp = "["+outp + "+"+str(startAtom.charge)+"]"

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
