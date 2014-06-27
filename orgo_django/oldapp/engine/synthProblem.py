from reactions import *
import randomGenerator
import random
import orgoStructure
import reactions as reactionsModule
import string
import serverRender
import copy
import itertools
try:
   import cPickle as pickle
except:
   import pickle


#ReactionStep class
    #Molecule box to point from
    #Molecule box to point to
    #List of added reagents
    #List of added molecules other than molecule box
class ReactionStep:
    def __init__(self, parentMoleculeBox):
        self.reactantBox = parentMoleculeBox
        self.otherMoleculeBox = MoleculeBox([])
        self.productBox = None
        
        #hasReagents is a dict that maps numbers (of reagents) to true/false
        self.hasReagents = parseReagentsString("")
    
    #Add more reagents.
    def addReagent(self, reagentDict):
        for reagent in list(REAGENTS):
            if self.hasReagents[reagent] or reagentDict[reagent]:
                self.hasReagents[reagent] = True
                
    #def addReagent(self, reagentBox): #reagentList is a list of reagent boxes
    #    for reagent in list(REAGENTS):
    #        if reagentBox.hasReagents[reagent]:
    #            self.hasReagents[reagent] = True
    #    return self.react()
          
    def addMolecule(self, moleculeBox):
        self.otherMoleculeBox = moleculeBox
        
    #Returns True or False depending on whether or not a reaction occurred.
    #If True, it updates self.product to be a new MoleculeBox containing 
    def react(self, mode="generate"):
        #mode can be "generate" or "check".  Generate returns true iff the reaction makes
        #*new* products; check returns true iff the reaction specified is a valid combination
        #of reagents.
        
        #a helper function
        def hasReagent(acceptableReagents):
            return True in [self.hasReagents[x] for x in acceptableReagents]
            
        for reaction in REACTIONS:
        #This loops through all of the reagents, checking if they are satisfied.
            if (True in [not hasReagent(acceptableReagents) for acceptableReagents in reaction[0]]): #returns False if the list comprehension [hasReagent...] is either empty or all full of Trues
                continue
            else:
                #react!
                try:
                    products = reaction[1](self.reactantBox.molecules)(self.otherMoleculeBox.molecules) #a function of two variables
                except ReactionTooCrazyError:
                    #TODO: write some sort of return that alerts the frontend.
                    return False
                if mode == "check":
                    #check if the output list is non-empty
                    #if so, reaction is successful
                    if products != []:
                        self.productBox = MoleculeBox(products)
                        return True
                    #if not, the old set of molecules remain intact
                    else:
                        self.productBox = MoleculeBox(self.reactantBox.molecules + self.otherMoleculeBox.molecules)
                        return True
                elif mode == "generate":
                    if len(products) > 4:
                        #Too many molecules.  It takes too long to moleculeCompare them pairwise, so just 
                        #return False now.
                        return False
                    #Return true if some new molecule was made during the course of the reaction
                    for product in products:
                        if sum([moleculeCompare(product, reactant) for reactant in self.reactantBox.molecules])==0:
                            self.productBox = MoleculeBox(products)
                            return True
                    return False
                else:
                    print "Invalid mode in react."
                    raise StandardError
                        
                
        return False
        
        #Check if reagent requirements are satisfied
        #Check if anything is produced by the reaction function
        
        
        
    #Returns a string of the reagents contained in this reaction step as proper HTML format.
    def stringList(self):
        out = ''
        for reagent in list(self.hasReagents):
            if self.hasReagents[reagent]:
                out += REAGENTS[reagent][0] + ', '
        return out
    
    
    def checkStep(self, target):
        #Takes in a molBox (target) and checks if its own products are identical to those in target.
        #Regardless of identical-ness, returns its own products as well.
        if self.react(mode="check"):
            return (boxEqualityChecker(self.productBox, target), self.productBox)
        else:
            return (boxEqualityChecker(self.reactantBox, target), self.reactantBox)
        #The input was ill-formatted.
        return (False, self.productBox)


#Called by checkIfEqualsTarget in MoleculeBoxModel in models
#Called by checkStep in ReactionStep in synthProblem
#Returns a boolean.
#first is a moleculebox and so is second
def boxEqualityChecker(first, second):
    assert isinstance(first, MoleculeBox)
    assert isinstance(second, MoleculeBox)
    #Does each product correspond to exactly one target?
    if len(first.molecules) != len(second.molecules):
        return False
    for output in first.molecules:
        OK = False
        for target in second.molecules:
            if moleculeCompare(output, target):
                OK = True
                second.molecules.remove(target)
                break
            #If by this point, we haven't found a match, return False.
        if not OK:
            return False
    #Reached the end of molecule list - must have perfect match.
    return True
    

#MoleculeBox class
#Represents a draggable box containing molecules (svg).
    #List of molecules contained in it
class MoleculeBox:
    def __init__(self, moleculesList):
        self.molecules = moleculesList

    #Returns a string of smiles of the molecules contained in this box, separated by spaces.
    def stringList(self):
        outp = ""
        for mol in self.molecules:
            outp += smiles(mol) + "."
        return serverRender.render(outp)
 
 

     
#Enter a string, such as "H2 cat Pd|C"
#Returns a dictionary, such as {"H2":True, "PDC":True, "ETOH":False, ...}
def parseReagentsString(inpstring):
    string = inpstring.lower()
    outp = {}
    
    if (string == ""):
        for reagent in list(REAGENTS):
            outp[reagent] = False
        return outp
        
    for reagent in list(REAGENTS):
        outp[reagent] = False
        for spelling in REAGENTS[reagent][1]:
            if spelling.lower() in string:
                outp[reagent] = True

    #no longer hacky
    #Check all valid reagent strings for whether they are substrings of OTHER valid reagent strings.
    #If so, check the count on the substring.
    #Make sure you don't count substrings if you're counting things they're part of.
    #We will precompute and load the precise set of things we need to check.
    
    #Every now and then, set recompute to True.
    recompute = True
    structure = []
    filename = "SubstringCheckerData.pickle"
    
    if recompute:
        def getOtherReagentStringsContaining(reagentID, substring):
            retval = []
            for key in REAGENTS:
                if not key==reagentID:
                    validStrings = REAGENTS[key][1]
                    for validString in validStrings:
                        validStringLower = validString.lower()
                        if substring in validStringLower:
                            retval += [validStringLower]
            return retval
        
        structure = []
        for key in REAGENTS:
            validStrings = REAGENTS[key][1]
            for validString in validStrings:
                validStringLower = validString.lower()
                get = getOtherReagentStringsContaining(key, validStringLower)
                if not (get == []):
                    structure += [(key, validStringLower, get)]
                    
                    
        ##write 'structure' down somewhere
        toWrite = pickle.dumps(structure)
        f = open(filename, 'w')
        f.write(toWrite)
        f.close()
                    
    else:
        ##read 'structure' from somewhere
        f = open(filename, 'r')
        toLoad = f.read()
        structure = pickle.loads(toLoad)
        f.close()
    
    for value in structure:
        val = value[0]
        sub = value[1]
        bigList = value[2]
        
        if (string.count(sub) != 0):
            if sum([string.count(big.lower()) for big in bigList]) >= sum([string.count(x.lower()) for x in REAGENTS[val][1]]):
                outp[val] = False
    
       
    return outp
    
def makeStartingMaterial(mode, count=1):
    #Make starting material, based on classes of reactions selected.
    molecules = []
    if ('10A Alkenes: halide addition' in mode) or ('10B Alkenes: other' in mode) or ('11 Alkynes' in mode):
        for i in xrange(count):
            if random.random() < 0.4:
                forceTerminalAlkyne = True
            else:
                forceTerminalAlkyne = False
            molecules.append(randomGenerator.randomStart(endProb=0.3, maxBranchLength=10,
            alkyneProb=0.1, alkeneProb=0.1,
            BrProb=0.1, ClProb=0.1, OHProb=0.05, forceTerminalAlkyne = forceTerminalAlkyne)[0])
        molecules = removeDuplicates(molecules)
        if debug:
            print "Starting material: " + str(smiles(molecules))
        return [MoleculeBox([molecule]) for molecule in molecules]
    
    
def randomSynthesisProblemMake(mode, steps = 20, maxLength = 30, count = 2):
    steps, fused = randomSynthesisProblemStart(mode, steps, maxLength, count)
    if fused:
        return steps
    productsNeeded = []  #Tracks all of the molecules we need.
    productsNeeded += steps[-1].reactantBox.molecules
    steps[-1].keep = True  #Add this new attribute to track which reactions are necessary.
    keepers = 0
    #If we found at least one more keeper in the last cycle, keep going.
    while sum([hasattr(step, "keep") for step in steps]) > keepers:
        keepers = sum([hasattr(step, "keep") for step in steps])
        for step in steps:
            for molecule in step.productBox.molecules:
                if molecule in productsNeeded:
                    step.keep = True
                    productsNeeded += (step.reactantBox.molecules + step.otherMoleculeBox.molecules)
    steps2 = copy.copy(steps)
    for step in steps:
        if not(hasattr(step, "keep")):
            steps2.remove(step)
    return steps2
    
    
def randomSynthesisProblemStart(mode, steps = 20, maxLength = 30, count=2):
    fused = False
    #Mode controls the reagents that are legal, as well as the distribution of starting materials.
    legalRxns = []
    for reactionSet in mode:
        legalRxns += typeToReaction[reactionSet]
    molBoxes = makeStartingMaterial(mode, count)
    reactions = [] #List of reactions that we want to keep
    
    #Each mode sets up its legal reactions and its starting materials.
    
    
    #Generate a list of legal addition (molecule-molecule bonding) reactions
    legalAddRxns = [reaction for reaction in legalRxns if 'add' in reaction[2]]
    #Give reactions marked "interesting" more emphasis, by cloning them.
    #Does anyone have a better idea for doing this?
    newLegalRxns = []
    for reaction in legalRxns:
        if 'interesting' in reaction[2]:
            newLegalRxns += [reaction]*5    #Change to increase/lessen emphasis
        else:
            newLegalRxns.append(reaction)
    legalRxns = newLegalRxns
 
    if debug:
        print "Initial mol:"
        for molBox in molBoxes:
            print smiles(molBox.molecules)
    
    #Try to react a bunch of times.
    for attemptNo in xrange(steps):
        #Tests for prematurely ending the generation process.
        if sum([len(molBox.molecules) for molBox in molBoxes]) > 4:
            if debug:
                print "Too many molecules!"
            return randomSynthesisProblemStart(mode, steps, maxLength, 1)
        for molBox in molBoxes:
            for molecule in molBox.molecules:
                if len(molecule.atoms) > maxLength:
                    if debug:
                        print "Molecule too large!"
                    return randomSynthesisProblemStart(mode, steps, maxLength, 1)
        newMolBoxes = []
        
        #Go through each molecule, and attempt a random reaction.
        for molBox in molBoxes:
            #There's a small chance of skipping. - Maybe delete, idk?
            if random.random() < .2:
                newMolBoxes.append(molBox)
                continue
            #Otherwise, come up with a random reaction, and try it out.
            reagents, rxnFunction, ignore = legalRxns[random.randint(0, len(legalRxns)-1)]
            if debug:
                print "Trying step: " + str(reagents)
            currentRxn = ReactionStep(molBox)
            for reagent in reagents:
                currentRxn.hasReagents[reagent[0]] = True
            if currentRxn.react() and len(currentRxn.productBox.molecules)+len(molBox.molecules)<5:
                #A good reaction.
                newMolBoxes.append(currentRxn.productBox)
                reactions.append(currentRxn)
                if debug:
                    print "Result of successful reaction: " +str(smiles(currentRxn.productBox.molecules))
 
            else:
                #Not a good reaction - that's OK, keep going
                newMolBoxes.append(molBox)
                
        molBoxes = []
        if len(newMolBoxes) == 1 or len(legalAddRxns)==0:
            #No point in trying to fuse molecules if you only have one molecule to begin with.
            molBoxes = newMolBoxes
            continue
        #Now, try to fuse molecules?
        reagents, rxnFunction, ignore = legalAddRxns[random.randint(0, len(legalAddRxns)-1)]
        if debug:
            print "Trying fusion: "+str(reagents)

        for i in xrange(len(newMolBoxes)):
            for j in xrange(i+1,len(newMolBoxes)):
                #Loop through all pairs of molecules.
                molBox1 = newMolBoxes[i]
                molBox2 = newMolBoxes[j]
                cc1 = molBox1.molecules[0].countElement('C')
                cc2 = molBox2.molecules[0].countElement('C')
                currentRxn = ReactionStep(molBox1)
                currentRxn.addMolecule(molBox2)
                for reagent in reagents:
                    currentRxn.hasReagents[reagent[0]] = True
                if currentRxn.react():
                    #OK, we have a reaction.  But, did we get fusion?
                    if sum([product.countElement('C') == cc1 + cc2 for product in currentRxn.productBox.molecules]) > 0:
                        #Success!
                        reactions.append(currentRxn)
                        newMolBoxes.remove(molBox1)
                        newMolBoxes.remove(molBox2)
                        molBoxes = newMolBoxes + [currentRxn.productBox]
                        fused = True
                        if debug:
                            print "Result: " +str(smiles(currentRxn.productBox.molecules))
                            print molBoxes
        if len(molBoxes) == 0:
            #Didn't fuse any molecules.  Oh well.
            molBoxes = newMolBoxes
    if len(reactions) == 0:
 
        return randomSynthesisProblemStart(mode, steps, maxLength, 1)
    return reactions, fused
                         
 
#def [moleculeboxes] = getStartingMoleculeBoxes(reactionSteps) in synthProblem
#Helper method used by a constructor in models.
def getStartingMoleculeBoxes(reactionSteps):
    products = list(set([reactionStep.productBox for reactionStep in reactionSteps]))
    allMolecules = list(set([step.reactantBox for step in reactionSteps]) | set([step.otherMoleculeBox for step in reactionSteps if (step.otherMoleculeBox.molecules != [])]))
    startingMoleculeBoxes = [molecule for molecule in allMolecules if (molecule not in products)]
    startingMoleculeBoxes = list(set(startingMoleculeBoxes)) #this should remove duplicates
    return startingMoleculeBoxes



def moleculeBoxHtml(moleculeBox):
    html = "<div class = \"molecule\" class=\"ui-widget-content\"  >"
    html += serverRender.render(moleculeBox.stringList())
    html += "</div>"
    return html   
 
def reactionStepHtml(reactionStep):
    html = ""
    for reagent in list(REAGENTS):
        if reactionStep.hasReagents[reagent]:
            html += REAGENTS[reagent][0] + ", "
            
    return "<div class = \"reaction\" class = \"ui-widget-content\">"+(html[:-2])+"<img src=\"http://felixsun.scripts.mit.edu/orgo/static/arrow.png\"/></div>"

                         
def generateNameReagentProblem(mode="AlkeneAlkyne"):
    #Endless loop, for now.  Maybe have some sort of give-up condition?
    while True:
        legalRxns = []
        for reactionSet in mode:
            legalRxns += typeToReaction[reactionSet]
        reactantBox = makeStartingMaterial(mode)[0]  #makeStartingMaterial returns a list of molBoxes
        #Try up to 10 times to make a reaction.
        for attemptNo in xrange(10):
            reagents, rxnFunction, labels = legalRxns[random.randint(0, len(legalRxns)-1)]
            if debug:
                print "Trying step: " + str(reagents)
            currentRxn = ReactionStep(reactantBox)
            for reagent in reagents:
                currentRxn.hasReagents[reagent[0]] = True
            if currentRxn.react() and len(currentRxn.productBox.molecules)<3:
                #A good reaction.
                #Decode labels to find the reaction label.
                for label in labels:
                    if not(label in nonLabelKeywords):
                        currentRxn.catagory = label
                return currentRxn
                if debug:
                    print "Result of successful reaction: " +str(smiles(currentRxn.productBox.molecules))
            else:
                #Try again.
                pass     





#Some constants
H2=1 
PDC=2
ETOH=3
HF=4
HBR=5
HCL=6
HI=7
CH2CL2=8
F2=9
BR2=10
CL2=11
I2=12
ROOR=13
RCO3H=14
H2SO4=15
H2O=16
HGSO4=17
BH3=18
THF=19
NAOH=20
H2O2=21
OSO4=22
NMO=23
ACETONE=24
O3=25
ME2S=26
ZN=27
LINDLAR=28
NA=29
NH3=30
NANH2=31
EQV1=32
HEAT=33
LIGHT=34
KOCCH33=35







#This gigantic terrible tuple is for determining which reaction should take place.

#First item of each tuple in this list:
    #a set of necessary reagents. Things listed together have an "or" relationship.
            #E.g. (("O3",),("CH2CL2",),("ME2S", "ZN")) means ozone AND ch2cl2 AND (me2s OR zn)
#Second item of each tuple in this list:
    #a function of two variables, which takes in a list of molecules (x) and another list of molecules (o) and returns them reacted
#These are listed roughly by precedence: earlier-listed reactions which qualify take precedence over later-listed ones.
REACTIONS = (
(((H2SO4,), (H2O,), (HGSO4,)), (lambda x: lambda o: acidhydrate(x+o, Molecule(Atom("O")), True)), ('10B Alkenes: other',)),
(((H2,),(PDC,),(ETOH,)), (lambda x: lambda o: hydrogenate(x+o)), ('10A Alkenes: halide addition',)),
(((HBR,), (ROOR,), (HEAT, LIGHT)), (lambda x: lambda o: radicalhydrohalogenate(x+o, "Br")), ('10A Alkenes: halide addition',)),
(((HBR,),(CH2CL2,), (EQV1,)), (lambda x: lambda o: hydrohalogenate1eq(x+o, "Br")), ('10A Alkenes: halide addition',)),
(((HF,),(CH2CL2,), (EQV1,)), (lambda x: lambda o: hydrohalogenate1eq(x+o, "F")), ('10A Alkenes: halide addition',)),
(((HI,),(CH2CL2,), (EQV1,)), (lambda x: lambda o: hydrohalogenate1eq(x+o, "I")), ('10A Alkenes: halide addition',)),
(((HCL,),(CH2CL2,), (EQV1,)), (lambda x: lambda o: hydrohalogenate1eq(x+o, "Cl")), ('10A Alkenes: halide addition',)),
(((HBR,),(CH2CL2,)), (lambda x: lambda o: hydrohalogenate(x+o, "Br")), ('10A Alkenes: halide addition',)),
(((HF,),(CH2CL2,)), (lambda x: lambda o: hydrohalogenate(x+o, "F")), ('10A Alkenes: halide addition',)),
(((HI,),(CH2CL2,)), (lambda x: lambda o: hydrohalogenate(x+o, "I")), ('10A Alkenes: halide addition',)),
(((HCL,),(CH2CL2,)), (lambda x: lambda o: hydrohalogenate(x+o, "Cl")), ('10A Alkenes: halide addition',)),
(((BR2,),(CH2CL2,),(EQV1,)), (lambda x: lambda o: halogenate1eq(x+o, "Br")), ('11 Alkynes',)),
(((F2,),(CH2CL2,),(EQV1,)), (lambda x: lambda o: halogenate1eq(x+o, "F")), ('11 Alkynes',)),
(((I2,),(CH2CL2,),(EQV1,)), (lambda x: lambda o: halogenate1eq(x+o, "I")), ('11 Alkynes',)),
(((CL2,),(CH2CL2,),(EQV1,)), (lambda x: lambda o: halogenate1eq(x+o, "Cl")), ('11 Alkynes',)),
(((BR2,),(CH2CL2,)), (lambda x: lambda o: halogenate(x+o, "Br")), ('10A Alkenes: halide addition',)),
(((F2,),(CH2CL2,)), (lambda x: lambda o: halogenate(x+o, "F")), ('10A Alkenes: halide addition',)),
(((I2,),(CH2CL2,)), (lambda x: lambda o: halogenate(x+o, "I")), ('10A Alkenes: halide addition',)),
(((CL2,),(CH2CL2,)), (lambda x: lambda o: halogenate(x+o, "Cl")), ('10A Alkenes: halide addition',)),
(((RCO3H,), (CH2CL2,)), (lambda x: lambda o: epoxidate(x+o)), ('10B Alkenes: other',)),
(((H2SO4,), (H2O,), (HGSO4,)), (lambda x: lambda o: acidhydrate(x+o, Molecule(Atom("O")), True)), ('10B Alkenes: other',)),
(((H2SO4,), (ETOH,), (HGSO4,)), (lambda x: lambda o: acidhydrate(x+o, ethanol, True)), ('10B Alkenes: other','illegal')),
(((H2SO4,), (HGSO4,)), (lambda x: lambda o: acidhydrate(x+o, x+o, True)), ('10B Alkenes: other','add')),
(((H2SO4,), (H2O,)), (lambda x: lambda o: acidhydrate(x+o, Molecule(Atom("O")))), ('10B Alkenes: other',)),
(((H2SO4,), (ETOH,)), (lambda x: lambda o: acidhydrate(x+o, ethanol)), ('10B Alkenes: other', 'illegal',)),
(((H2SO4,),), (lambda x: lambda o: acidhydrate(x+o, x+o)), ('10B Alkenes: other',)),
(((BR2,), (H2O,)), (lambda x: lambda o: halohydrate(x+o, Molecule(Atom("O")), "Br")),('10A Alkenes: halide addition',)),
(((BR2,), (ETOH,)), (lambda x: lambda o: halohydrate(x+o, ethanol, "Br")), ('10B Alkenes: other','illegal',)),
(((BR2,),), (lambda x: lambda o: halohydrate(x+o, x+o, "Br")), ('10B Alkenes: other',)),
(((I2,), (H2O,)), (lambda x: lambda o: halohydrate(x+o, Molecule(Atom("O")), "I")),('10A Alkenes: halide addition',)),
(((I2,), (ETOH,)), (lambda x: lambda o: halohydrate(x+o, ethanol, "I")), ('10B Alkenes: other','illegal')),
(((I2,),), (lambda x: lambda o: halohydrate(x+o, x+o, "I")), ('10B Alkenes: other',)),
(((F2,), (H2O,)), (lambda x: lambda o: halohydrate(x+o, Molecule(Atom("O")), "F")), ('10A Alkenes: halide addition',)),
(((F2,), (ETOH,)), (lambda x: lambda o: halohydrate(x+o, ethanol, "F")), ('10B Alkenes: other','illegal')),
(((F2,),), (lambda x: lambda o: halohydrate(x+o, x+o, "F")), ('10B Alkenes: other',)),
(((CL2,), (H2O,)), (lambda x: lambda o: halohydrate(x+o, Molecule(Atom("O")), "Cl")),('10A Alkenes: halide addition',)),
(((CL2,), (ETOH,)), (lambda x: lambda o: halohydrate(x+o, ethanol, "Cl")), ('10B Alkenes: other','illegal')),
(((CL2,),), (lambda x: lambda o: halohydrate(x+o, x+o, "Cl")), ('10B Alkenes: other',)),
(((BH3,), (THF,), (NAOH,), (H2O2,)), (lambda x: lambda o: hydroborate(x+o)), ('10B Alkenes: other',)),
(((BH3,), (THF,)), (lambda x: lambda o: hydroborate1(x+o)), ('10B Alkenes: other',)),
(((NAOH,), (H2O2,)), (lambda x: lambda o: hydroborate2(x+o)), ('10B Alkenes: other',)),
(((OSO4,), (NMO,), (ACETONE, H2O)), (lambda x: lambda o: dihydroxylate(x+o)), ('10B Alkenes: other',)),
(((O3,),(CH2CL2,),(ME2S,ZN)), (lambda x: lambda o: ozonolyse(x+o)), ('10B Alkenes: other',)),
(((NA,), (NH3,)), (lambda x: lambda o: sodiumAmmonia(x+o)), ('11 Alkynes',)),
(((LINDLAR,), (H2,)), (lambda x: lambda o: lindlar(x+o)), ('11 Alkynes',)),
(((NANH2,), (NH3,)), (lambda x: lambda o: alkyneDeprotonate(x+o)), ('11 Alkynes',)),
(((KOCCH33,),), (lambda x: lambda o: tertButoxide(x+o)), ('10B Alkenes: other', 'interesting')),
((), (lambda x: lambda o: acetylideAdd(x+o, x+o)),('11 Alkynes','add'))
)


#SYNTHONLY = [
#(((NANH2,), (NH3,)), (lambda x: lambda o: acetylideAdd(alkyneDeprotonate(x), o)),('11 Alkynes','add')),
#(((NANH2,), (NH3,)), (lambda x: lambda o: acetylideAdd(x, alkyneDeprotonate(o))),('11 Alkynes','add'))
#]
SYNTHONLY = [
(((NANH2,), (NH3,)), (lambda x: lambda o: acetylideAdd(alkyneDeprotonate(x+o), x+o)),('11 Alkynes','add')),
(((NANH2,), (NH3,)), (lambda x: lambda o: acetylideAdd(x+o, alkyneDeprotonate(x+o))),('11 Alkynes','add'))
]

FORSYNTH = [reaction for reaction in REACTIONS if not ('illegal' in reaction[2])] + SYNTHONLY
nonLabelKeywords = ['illegal', 'add', 'interesting']
#Make a dictionary mapping each reaction catagory to all reactions in that catagory.
typeToReaction = dict()
for reaction in FORSYNTH:
    for label in reaction[2]:
        if label in nonLabelKeywords:
            continue
        if label in typeToReaction:
            typeToReaction[label].append(reaction)
        else:
            typeToReaction[label] = [reaction]

#You can make larger groups of reactions down here.
AlkeneAlkyneMode = ('10A Alkenes: halide addition', '10B Alkenes: other', '11 Alkynes')



#This dictionary is for understanding what people type in.

#List of synonyms of reagents.
#The first value in the tuple is how it should be printed (with proper subscripts).
#The second value in the tuple is all of the items that should be recognizable as user input.
#The third value in the tuple is any of the items that the dropdown menu wants to display.
    #Preferably, make sure that the HTML output (simplified) is always included in the dropdown menu tuple.

#When adding a new reagent, also be sure to update REACTIONS, to add a value for the constant,
#and to update the section flagged "hacky" in the methods above.
#Also update the drop-down field in the frontend with the new typeable values.
REAGENTS = {
H2: ("H<sub>2</sub>",("H2", "Hydrogen"), ("H2",)),
PDC: ("Pd|C", ("Pd/C", "PdC", "Pd|C", "Pd C", "Palladium"), ("Pd/C",)),
ETOH: ("EtOH", ("EtOH", "Ethanol", "Ethyl alcohol", "C2H5OH"), ("EtOH",)),
HF: ("HF", ("HF", "Hydrogen fluoride", "Hydrofluoric acid"), ("HF",)),
HBR: ("HBr", ("HBr", "Hydrogen bromide", "Hydrobromic acid"), ("HBr",)),
HCL: ("HCl", ("HCl", "Hydrogen chloride", "Hydrochloric acid"), ("HCl",)),
HI: ("HI", ("HI", "Hydrogen iodide", "Hydroiodic acid"), ("HI",)),
CH2CL2: ("CH<sub>2</sub>Cl<sub>2</sub>", ("CH2Cl2", "Dichloromethane"), ("CH2Cl2",)),
F2: ("F<sub>2</sub>", ("Fluorine", "F2"), ("F2",)),
BR2: ("Br<sub>2</sub>", ("Bromine", "Br2"), ("Br2",)),
CL2: ("Cl<sub>2</sub>", ("Chlorine", "Cl2"), ("Cl2",)),
I2: ("I<sub>2</sub>", ("Iodine", "I2"), ("I2",)),
ROOR: ("ROOR", ("ROOR", "tBuOOtBu", "Peroxide"), ("ROOR","Peroxide")),
RCO3H: ("RCO<sub>3</sub>H",("mCPBA", "PhCO3H", "RCO3H"), ("mCPBA","PhCO3H","RCO3H")),
H2SO4: ("H<sub>2</sub>SO<sub>4</sub>", ("H2SO4", "Sulfuric acid", "H+"), ("H2SO4",)),
H2O: ("H<sub>2</sub>O", ("H2O", "Water", "HOH", "dihydrogen monoxide"), ("H2O",)),
HGSO4: ("HgSO<sub>4</sub> accels.", ("HgSO4", "Hg2+", "Mercury sulfate", "Mercury(II", "Mercury (II", "Mercury(2", "Mercury (2"), ("HgSO4",)),
BH3: ("BH<sub>3</sub>", ("BH3", "Borane"), ("BH3",)),
THF: ("THF", ("THF", "Tetrahydrofuran"), ("THF",)),
NAOH: ("NaOH", ("NaOH", "Sodium hydroxide", "Hydroxide", "OH-"), ("NaOH",)),
H2O2: ("H<sub>2</sub>O<sub>2</sub>", ("H2O2", "Hydrogen peroxide"), ("H2O2",)),
OSO4: ("OsO<sub>4</sub>", ("OsO4", "Osmium tetroxide", "Osmium oxide"), ("OsO4",)),
NMO: ("NMO", ("NMO", "NMMO", "N-Methylmorpholine N-oxide"), ("NMO",)),
ACETONE: ("Acetone", ("Acetone", "Propanone", "(CH3)2CO"), ("Acetone",)),
O3: ("O<sub>3</sub>", ("Ozone", "O3"), ("O3",)),
ME2S: ("Me<sub>2</sub>S", ("Dimethyl sulfide", "Methylthiomethane", "Me2S"), ("Me2S",)),
ZN: ("Zn", ("Zn", "Zinc"), ("Zn",)),
LINDLAR: ("cat. Lindlar", ("Lindlar",), ("Lindlar",)),
NA: ("Na", ("Sodium", "Na"), ("Na",)),
NH3: ("NH<sub>3 (L)</sub>", ("NH3", "Ammonia"), ("NH3",)),
NANH2: ("NaNH<sub>2</sub>", ("Amide", "Sodium amide", "NaNH2"), ("NaNH2",)),
EQV1: ("1 equiv.", ("1", "one"), ("1 equiv.",)),
HEAT: ("Heat", ("heat", "delta", "hot", "warm"), ("Heat",)),
LIGHT: ("Light", ("hv", "light", "hnu", "h v", "h nu", "photon"), ("Light",)),
KOCCH33: ("Tert-butoxide", ("Tert-butoxide", "OC(CH3)3", "OtBu", "tertbutoxide", "tert butoxide"), ("Tert-butoxide",))
}


#Debugging
if __name__ == "__main__":
    for reaction in REACTIONS:
        out = ""
        for reagent in reaction[0]:
            out += REAGENTS[reagent[0]][1][0] + " "
        print out



def reagentAutocompleteMake():
    historicalValue = '''['H2', 'Hydrogen', 'Pd/C', 'Palladium/Carbon catalyst', 'EtOH', 'Ethanol', 'Ethyl alcohol', 'C2H5OH', 'HF', 'Hydrogen fluoride', 'Hydrofluoric acid', 'HBr', 'Hydrogen bromide', 'Hydrobromic acid', 'HCl', 'Hydrogen chloride', 'Hydrochloric acid', 'HI', 'Hydrogen iodide', 'Hydroiodic acid', 'CH2Cl2', 'Dichloromethane', 'Fluorine', 'F2', 'Bromine', 'Br2', 'Chlorine', 'Cl2', 'Iodine', 'I2', 'ROOR', 'tBuOOtBu', 'Peroxide', 'Tert-butyl peroxide', 'Di-tert-butyl peroxide', 'mCPBA', 'PhCO3H', 'RCO3H', 'H2SO4', 'Sulfuric acid', 'H2O', 'Water', 'HOH', 'H20', 'HgSO4', 'Hg2+', 'Mercury sulfate', 'BH3', 'Borane', 'THF', 'Tetrahydrofuran', 'NaOH', 'Sodium hydroxide', 'Hydroxide', 'OH-', 'H2O2', 'Hydrogen peroxide', 'OsO4', 'Osmium tetroxide', 'Osmium oxide', 'NMO', 'NMMO', 'N-Methylmorpholine N-oxide', 'Acetone', 'Propanone', '(CH3)2CO', 'Ozone', 'O3', 'Dimethyl sulfide', 'Methylthiomethane', 'Me2S', 'Zn', 'Zinc', 'Lindlar catalyst', 'cat. Lindlar', 'Sodium', 'Na', 'NH3', 'Ammonia', 'Sodium amide', 'Sodamide', 'NaNH2', 'Amide', '1 equivalent', 'One equivalent', 'heat', 'hv', 'light', 'hnu', 'tert-butoxide', 'KOtBu', 'Potassium tert-butoxide', 'KOC(CH3)3']''' 
    #return historicalValue
    
    output = []
    
    for reagent in REAGENTS.itervalues():
        for typable in reagent[2]:
            output += [typable]
    
    return output
    
def reactionAutocompleteMake():

    ##Make this True every now and then for brief periods of time.
    recompute = True
    filename = "reactionAutocomplete.txt"
    
    ##reads precomputed constant thing from a file
    if not recompute:

        f = open(filename, 'r')
        output = f.read()
        f.close()
        
        return output
    
    ###Try to run the below computation every now and then, thereby replacing the output to historicalValue.
    if recompute:
        verbose = False
        listOfReactionsByReagent = []
        
        for reaction in REACTIONS:
            if verbose:
                print "reaction="+str(reaction)
                print "listOfReactionsByReagent"+str(listOfReactionsByReagent)
            reagentsToAdd = []
            reagents = reaction[0] #e.g. reagents = ((HBR,), (ROOR,), (HEAT, LIGHT))
            #We want all possible n-length reagent permutations, so that people can type in any order.
            n = len(reagents)
            for indices in [range(n)]: #itertools.permutations(range(n), n): # #e.g. indices = (0,1,2)
                if verbose:
                    print "indices="+str(indices)
                    print "reagentsToAdd="+str(reagentsToAdd)
                listInProgress = [""]     #is a list of strings
                for index in indices:   #e.g. index=1
                    if verbose:
                        print "index="+str(index)
                        print "listInProgress="+str(listInProgress)
                    tempListInProgress = listInProgress
                    listInProgress = []
                    for item in tempListInProgress: #item is a string
                        if verbose:
                            print "item="+str(item)
                        for reagent in reagents[index]: #e.g. reagent=LIGHT
                            if verbose:
                                print "reagent="+str(reagent)
                            for typable in REAGENTS[reagent][2]: #e.g. REAGENTS[LIGHT][2] = ("Light",)
                                if verbose:
                                    print "typable="+str(typable)
                                if item == "":
                                    listInProgress += [typable]
                                else:
                                    listInProgress += [(item + ", " + typable)]
                reagentsToAdd += listInProgress
            listOfReactionsByReagent += reagentsToAdd
            
        
        output =  str(listOfReactionsByReagent)
        
        #write it to a file
        f = open(filename, 'w')
        f.write(output)
        f.close()
        
        #return it so the website doesn't break
        return output



