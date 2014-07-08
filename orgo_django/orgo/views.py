# Create your views here.
from django.shortcuts import render
from django.contrib.auth import *
from django.contrib.auth.models import User
from django.contrib.auth.forms import PasswordResetForm, PasswordChangeForm
from django.contrib.auth.decorators import login_required
from django.http import HttpResponse
from django.views.decorators.csrf import csrf_exempt
from django.core.mail import EmailMessage
from django.utils import timezone
from django.utils.html import escape
import json
import traceback
import forkit
import engine.renderSVG as serverRender
import engine.reactions as reactions
import engine.randomGenerator as randomGenerator
import engine.molecularStructure as orgoStructure
from engine.toSmiles import smilesify
import orgo.models as models
from engine.synthProblem import *
TIMEOUT = 10 #seconds

def home(request, debug = ""):
    #Home page.
    if request.user.is_authenticated():
        return loggedInHome(request)
    outSmiles = smilesify(randomGenerator.randomStart()[0])
    svg = serverRender.render(outSmiles)
    return render(request, 'index.html', {'molecule': svg,
                                          'signUpForm': models.mySignUpForm,
                                          'logInForm': models.myAuthenticationForm(),
                                          'debug':debug , 
                                          'resetPWForm': models.myPasswordResetForm()})

def signUp(request):
    if request.method == 'POST':
        form = models.mySignUpForm(request.POST)
        if form.is_valid():
            if User.objects.filter(email=form.cleaned_data['email']):
                #Non-unique email - not allowed!
                return home(request, debug="Someone else has already registered with this email!")
            form.save()
            newUser = authenticate(username = form.cleaned_data['username'], 
                password = form.cleaned_data['password1'])
            login(request, newUser)
            return loggedInHome(request)
        else:
            return home(request, debug = "<span style=\"color:FF0000\">Bad account info, please try again.</span>")
            
def logIn(request):
    if request.method == 'POST':
        username = request.POST['username']
        password = request.POST['password']
        #Do log in check stuff
        user = authenticate(username = username, password = password)
        if user != None:
            login(request, user)
            return loggedInHome(request)
    return home(request, debug = "<span style=\"color:FF0000\">Invalid login, sorry.</span>")
    


def logOut(request):
    logout(request)
    return home(request)

def resetPW(request):
    if request.method == 'POST':
        #Use built-in Django magic to reset a user's password and send him an email.
        email = request.POST['email']
        try:
            user = User.objects.get(email=email)
        except:
            return home(request, debug = "The email you entered does not correspond to any user.")
        password = User.objects.make_random_password()
        user.set_password(password)
        user.save()
        subject = "Orgo - your new password"
        body = 'Hello '+user.username+''',
Here is your new password: ''' + password + '''
Because this password was sent over email, it is not secure as a permanent password.  Please log in and change it immediately.
Thank you,
Your friendly admins'''
        newMessage = EmailMessage(subject, body, to=[email])
        newMessage.send()
        return render(request, 'successfulReset.html')
    return home(request, debug = "Unknown error in password reset")
    
    
def changePW(request):
    #Change the user's password.
    if request.method == 'POST':
        if authenticate(username=request.user.username, password=request.POST['old_password']) == None:
            return loggedInHome(request, debug="Your old password is incorrect.")
        if request.POST['new_password1'] != request.POST['new_password2']:
            return loggedInHome(request, debug="Your two new passwords don't match.")
        request.user.set_password(request.POST['new_password1'])
        request.user.save()
        return loggedInHome(request, debug="Password changed successfully.")
        #old_password, new_password1, new_password2


#@login_required
#def returnToLoggedInHome(request):
#    return loggedInHome(request)

    

###Can delete; this is me learning Django
def outpSmiles(request):
    if request.method == 'POST':
        form = models.MoleculeForm(request.POST)
        if form.is_valid():
            molecule = models.MoleculeModel.create(request.POST['smiles'])
            molecule.save()
            return renderSmiles(request, molecule)
    
    
@csrf_exempt
def homeMoleculeChanger(request):
    #Returns new molecules for the AJAX tester (home molecule changer)
    outSmiles = smilesify(randomGenerator.randomStart()[0])
    svg = serverRender.render(outSmiles)
    return HttpResponse(svg)
            
@login_required
def loggedInHome(request, debug = ""):
    #Home page for those who have logged in.
    profile = request.user.profile
    #Set up the ChooseReagentsForm with the last choices the user made.
    initialValuesDict = dict()
    graphList = []
    for name, reactions in typeToReaction.items():
        try:
            profile.savedReagentTypes.get(name=name)
        except:
            initialValuesDict[name] = False
        else:
            initialValuesDict[name] = True
        #Load up the graph data
        try:
            accuracyObj = profile.accuracies.get(catagory=name)
        except:
            pass
        else:
            graphList.append([1.0*accuracyObj.correct/accuracyObj.total, name, accuracyObj.correct, accuracyObj.total])
    graphList.sort(key=lambda item: item[1])
    
    #Load high scores.
    bestUsers = models.UserProfile.objects.order_by('-correctSynths', 'user__username')[:5]
    synthHighScore = ""
    for i in xrange(min(5, len(bestUsers))):
        prof = bestUsers[i]
        if prof==profile:
            synthHighScore += "<tr><td>"+str(i+1)+"</td><td><b>"+prof.user.username+"</b></td><td>"+str(prof.correctSynths)+"</td></tr>"
        else:
            synthHighScore += "<tr><td>"+str(i+1)+"</td><td>"+prof.user.username+"</td><td>"+str(prof.correctSynths)+"</td></tr>"
    if request.user.profile not in bestUsers:
        #Insert a gap in the table
        synthHighScore += "<tr><td></td><td>...</td><td></td></tr>"
        #Find rank of user.
        for index, prof in enumerate(models.UserProfile.objects.order_by('-correctSynths', 'user__username')):
            if prof == profile:
                synthHighScore += "<tr><td>"+str(index+1)+"</td><td><b>"+request.user.username+"</b></td><td>"+str(profile.correctSynths)+"</td></tr>"
        synthHighScore += "<tr><td></td><td>...</td><td></td></tr>"
        
    #Load most helpful.
    mostHelpful = models.UserProfile.objects.order_by('-assists', 'user__username')[:5]
    helpHighScore = ""
    for i in xrange(min(5, len(mostHelpful))):
        prof = mostHelpful[i]
        if prof==profile:
            helpHighScore += "<tr><td>"+str(i+1)+"</td><td><b>"+prof.user.username+"</b></td><td>"+str(prof.assists)+"</td></tr>"
        else:
            helpHighScore += "<tr><td>"+str(i+1)+"</td><td>"+prof.user.username+"</td><td>"+str(prof.assists)+"</td></tr>"
    if request.user.profile not in mostHelpful:
        #Insert a gap in the table
        helpHighScore += "<tr><td></td><td>...</td><td></td></tr>"
        #Find rank of user.
        for index, prof in enumerate(models.UserProfile.objects.order_by('-assists', 'user__username')):
            if prof == profile:
                helpHighScore += "<tr><td>"+str(index+1)+"</td><td><b>"+request.user.username+"</b></td><td>"+str(profile.assists)+"</td></tr>"
        helpHighScore += "<tr><td></td><td>...</td><td></td></tr>"
    return render(request, 'loggedin.html', {'name': request.user.username, 
                                             'ChooseReagentsForm':models.ChooseReagentsForm(initial=initialValuesDict), 
                                             'debug': debug,
                                             'graphData': graphList,
                                             'changePW': PasswordChangeForm(request.user),
                                             'synthHighScore': synthHighScore,
                                             'helpHighScore': helpHighScore})
    
###Can delete; this is me learning Django
def renderSmiles(request, molecule):
    return render(request, 'loggedin.html', {'item1': serverRender.render(str(molecule.smiles))})

    


##Figure out how to generate all possible permutations of allowable reactions for this.    
reagentAutocomplete=reagentAutocompleteMake()  #defined in synthProblem
reactionAutocomplete=reactionAutocompleteMake()  #defined in synthProblem


    
@login_required    
def renderOldNameReagent(request):
    profile = request.user.profile
    try:
        step = profile.currentNameReagentProblem
    except:
        return renderNameReagent(request)
    autocompleteType = profile.autocompleteType
    if autocompleteType == "Reactions":
        autocomplete = reactionAutocomplete
    elif autocompleteType == "Reagents":
        autocomplete = reagentAutocomplete
    elif autocompleteType == "None":
        autocomplete = ""
    else:
        return loggedInHome(request, debug = "An invalid autocomplete mode was picked.")
    return render(request, 'problemInterface.html', {"ReactantMolecule": step.reactantBox.svg, 
                                                     "TargetMolecule": step.productBox.svg, 
                                                     "Name": request.user.username,
                                                     "Autocomplete": autocomplete,
                                                     "needsHelp": 0,})
    

def renderNameReagent(request, tutorial=False):
    try:
        profile = request.user.profile
    except:
        #Anonymous user.
        problem = generateNameReagentProblem(AlkeneAlkyneMode)
        step = models.ReactionStepModel.create(problem)
        step.save()
        request.session['problem'] = step
        autocomplete = reactionAutocomplete
        if 'tutorial' in request.session:
            tutorial = False
        else:
            tutorial = True
            request.session['tutorial'] = True
    else:
        #Sometimes, the user doesn't even have a previous problem, so deleting doesn't always work.
        try:
            temp1 = profile.currentNameReagentProblem.productBox
            temp2 = profile.currentNameReagentProblem.reactantBox
            profile.currentNameReagentProblem.delete()
            temp1.delete()
            temp2.delete()
        except:
            pass
            
           
        modes, autocomplete = checkboxUpdate(request, profile)
        if modes == []:
            #Error - at least one mode must be selected!
            return loggedInHome(request, debug = "You must pick at least one reaction type!")
        
        problem = generateNameReagentProblem(modes)
        step = models.ReactionStepModel.create(problem)
        step.save()
        profile.currentNameReagentProblem = step
        profile.save()
    return render(request, 'problemInterface.html', {"ReactantMolecule": step.reactantBox.svg,
                                                     "TargetMolecule": step.productBox.svg, 
                                                     "Name": request.user.username,
                                                     "Autocomplete": autocomplete,
                                                     "needsHelp": int(tutorial),})
        
        
def checkboxUpdate(request, profile): 
    #profile = request.user.profile
    modes = []
    if request.method == 'POST':
        #User filled out the checkboxes
        #First, clear the user's savedReagentTypes
        profile.savedReagentTypes.clear()
        #Parse and save which reaction types were checked.
        checkboxes = models.ChooseReagentsForm(request.POST)
        if checkboxes.is_valid():
            for name, reactions in typeToReaction.items():
                if checkboxes.cleaned_data[name]:
                    modes.append(name)
                    try:
                        reagentType = models.ReagentType.objects.get(name = name)
                    except:
                        #Make a new reagentType.  (We only need to do this once after we update our reagent types.)
                        reagentType = models.ReagentType.create(name = name)
                        reagentType.save()
                    profile.savedReagentTypes.add(reagentType)
       
            autocompleteType = checkboxes.cleaned_data["autocomplete"]
            if autocompleteType == "Reactions":
                autocomplete = reactionAutocomplete
            elif autocompleteType == "Reagents":
                autocomplete = reagentAutocomplete
            elif autocompleteType == "None":
                autocomplete = ""
            else:
                return loggedInHome(request, debug = "An invalid autocomplete mode was picked.")
            profile.autocompleteType = autocompleteType
            profile.save()
        else:
            return ([], None)
    else:
        #User clicked on "new problem"
        #Load up the user's reagent preferences.
        for reagentTypeInstance in profile.savedReagentTypes.all():
            modes.append(reagentTypeInstance.name)
        autocompleteType = profile.autocompleteType
        if autocompleteType == "Reactions":
            autocomplete = reactionAutocomplete
        elif autocompleteType == "Reagents":
            autocomplete = reagentAutocomplete
        elif autocompleteType == "None":
            autocomplete = ""
        else:
            return ([], None)
    return modes, autocomplete


    
@csrf_exempt
def checkNameReagent(request):
    #Takes in a list of reagents that the user guessed.
    #Returns whether that list is correct, plus the actual product, if that list is incorrect.
    if request.method == 'POST':
      try:
        reagentsDict = parseReagentsString(request.POST['reagents'])
        try:
            problem = request.user.profile.currentNameReagentProblem
            anonymous = False
        except:
            #User is anonymous.
            problem = request.session['problem']
            anonymous = True
        reactant = problem.reactantBox.moleculeBox
        target = problem.productBox.moleculeBox
        testStep = ReactionStep(reactant)
        #If the problem is done, do nothing - don't let them get more points.
        if problem.done:
            return HttpResponse("")
        testStep.addReagent(reagentsDict)
        correct, products = testStep.checkStep(target)
        responseData = dict()
        responseData["success"] = correct
        if isinstance(products, str):
            responseData["product"] = products
        else:
            responseData["product"] = products.stringList()
        if not anonymous:
            profile = request.user.profile
            #If we have the correct answer, free up some database space by deleting this stuff.
            thisCatagory = profile.currentNameReagentProblem.catagory
            try:
                accModel = profile.accuracies.get(catagory=thisCatagory)
            except:
                accModel = models.AccuracyModel.create(catagory=thisCatagory)
                accModel.save()
                profile.accuracies.add(accModel)
            if correct == True:
                accModel.correct += 1
                accModel.total += 1
                accModel.save()
                problem.done = True
                problem.save()
            else:
                accModel.total += 1
                accModel.save()
        else:
            if correct == True:
                problem.done = True
                problem.save()
      except StandardError as e:
        responseData = dict()
        responseData["product"] = str(e)
        responseData["success"] = False
        
      return HttpResponse(json.dumps(responseData))
      
@csrf_exempt
def showNRAnswer(request):
    try:
        problem = request.user.profile.currentNameReagentProblem
    except:
        #Anonymous user
        problem = request.session['problem']
    out = problem.html
    htmlOutput =  "<li class=\"reagent\" class = \"ui-state-default\" >"+out+"<img src=\"http://felixsun.scripts.mit.edu/orgo/static/arrow.png\"/></li>"
    #Mark this problem as done.
    problem.done = True
    problem.save()
    return HttpResponse(htmlOutput)
 
    
##Make this have a shiny flowchart layout once that becomes possible.
##Obsolete? Not convinced this is actually used anywhere.
# def moleculesAndReactionsHtml(startingMaterials, reactionSteps):
    # html = ""
    # startingMaterialsCopy = copy.deepcopy(startingMaterials)
    # for startingMaterial in startingMaterialsCopy:
        # keepgoing = True
        # moleculeToCheck = startingMaterial
        # while keepgoing:
            # html += moleculeBoxHtml(moleculeToCheck)
            # stepList = [reactionStep for reactionStep in reactionSteps if reactionStep.reactantBox == moleculeToCheck]
            # if len(stepList) == 1:
                # html += reactionStepHtml(stepList[0])
            # elif len(stepList) == 0:
                # keepgoing = False
                # html += "<br/>"
            # else:
                # html += reactionStepHtml(stepList[0])
                
    
    # return html
    
@csrf_exempt
def makeReagentHtml(request):
    try:
        reagentString = request.POST['reagentString']
        dict = parseReagentsString(reagentString)
        html = ""
        for reagent in dict:
            if dict[reagent]:
                html += REAGENTS[reagent][0] + ", "
        if html == "":
            htmlOutput = "<li class=\"reagent\" class = \"ui-state-default\" reagentString = \""+reagentString+"\">[no valid reagents]<img src=\"http://felixsun.scripts.mit.edu/orgo/static/arrow.png\"/></li>"
        else:
            htmlOutput =  "<li class=\"reagent\" class = \"ui-state-default\" reagentString = \""+reagentString+"\">"+html[:-2]+"<img src=\"http://felixsun.scripts.mit.edu/orgo/static/arrow.png\"/></li>"
        return HttpResponse(htmlOutput)
    except:
        tb = traceback.format_exc()
        return HttpResponse(str(tb))
        
###Some methods for displaying and operating with giant synthesis problems.
@login_required
def renderOldSynthesis(request):
    profile = request.user.profile
    try:
        synthesis = profile.currentSynthesisProblem
    except:
        return renderSynthesis(request)
    autocompleteType = profile.autocompleteType
    if autocompleteType == "Reactions":
        autocomplete = reactionAutocomplete
    elif autocompleteType == "Reagents":
        autocomplete = reagentAutocomplete
    elif autocompleteType == "None":
        autocomplete = ""
    else:
        return loggedInHome(request, debug = "An invalid autocomplete mode was picked.")
    return render(request, 'synthesisProblemInterface.html', {"TargetMolecule": synthesis.target.svg, 
                                                              "Name": request.user.username,
                                                              "Autocomplete": autocomplete,
                                                              "needsHelp": 0})

def loadSynthesisFromId(request):
    #Make a deep copy of the problem, and render it to the new user.
    if request.method != 'POST':
        #This should never happen.
        return
    loadId = request.POST['Id']
    try:
        synthesis = models.SynthesisProblemModel.objects.get(pk=loadId)
    except:
        return loggedInHome(request, debug="Invalid ID number.")

    #We may need to unassociate this problem from the old user.
    tempManager = synthesis.userprofile_set.all()
    if len(tempManager) != 0:
        oldUserP = tempManager[0]
        oldUserP.currentSynthesisProblem = None
        oldUserP.save()
    sCopy = forkit.tools.fork(synthesis, deep=True)
    sCopy.retain = False
    sCopy.save()
    profile = request.user.profile
    profile.currentSynthesisProblem = sCopy
    profile.save()
    return renderOldSynthesis(request)

@login_required
def renderSynthesis(request, tutorial = False):
    profile = request.user.profile
    try:
        profile.currentSynthesisProblem
    except:
        #Oops, accidentally deleted it.
        pass
    else:
        #If the retain attribute is True, don't delete.
        if profile.currentSynthesisProblem == None or not(profile.currentSynthesisProblem.retain):
            #Sometimes, the user doesn't even have a previous problem, so deleting doesn't always work.
            try:
                for arrowModel in profile.currentSynthesisProblem.arrows.all():
                    arrowModel.delete()
                for moleculeModel in profile.currentSynthesisProblem.molecules.all():
                    moleculeModel.delete()
                for arrowModel in profile.currentSynthesisProblem.solution.arrows.all():
                    arrowModel.delete()
                for moleculeModel in profile.currentSynthesisProblem.solution.molecules.all():
                    moleculeModel.delete()
                profile.currentSynthesisProblem.solution.delete()
                profile.currentSynthesisProblem.target.delete()
                profile.currentSynthesisProblem.delete()
            except:
                pass
    
    modes, autocomplete = checkboxUpdate(request, profile)
    if modes == []:
        #Error - at least one mode must be selected!
        return loggedInHome(request, debug = "You must pick at least one reaction type!")
    
    assert modes != []
    reactionsteps = randomSynthesisProblemMake(modes)
    synthesis = models.SynthesisProblemModel.create(reactionsteps)
    synthesis.save()
    profile.currentSynthesisProblem = synthesis
    profile.save()
    return render(request, 'synthesisProblemInterface.html', {"TargetMolecule": synthesis.target.svg, 
                                                              "Name": request.user.username,
                                                              "Autocomplete": autocomplete,
                                                              "needsHelp": int(tutorial)})

    

def getSynthesisData(request, synthesis=None):
    #Should return a JSON string with the following attributes contained:
    #molecules is an array of arrays: [ [idnumber, "<svg>...</svg>"], ... ]
    #arrows is an array of arrays: [ [idnumber1, idnumber2, "reagentText"], ...]
    #success   -- a boolean (true/false)
    
    #The "Return dictionary" mode is for the helper to see the same problem.
    if synthesis == None:
        mode = "Send to browser"
        synthesis = request.user.profile.currentSynthesisProblem
    else:
        mode = "Return dictionary"
    
    #Iterate over all molecules for a specific synthesis
    try:
        moleculesOutput = [ (moleculeBoxModel.id, moleculeBoxModel.svg) 
                        for moleculeBoxModel in synthesis.molecules.all()]
    
        #Iterate over all arrows for a specific synthesis
        
        arrowsOutput = []
        for arrowModel in synthesis.arrows.all():
            assert arrowModel.pointFrom != None
            assert arrowModel.pointTo != None
            arrowsOutput += [(arrowModel.pointFrom.id, arrowModel.pointTo.id, arrowModel.reagentsHtml)]
            
        arrowsOutput = [ (arrowModel.pointFrom.id, arrowModel.pointTo.id, arrowModel.reagentsHtml)
                         for arrowModel in synthesis.arrows.all()]
        
        responseData = dict()
        
        responseData["success"] = synthesis.checkIfSolved()
        responseData["molecules"] = moleculesOutput
        responseData["arrows"] = arrowsOutput
        if mode == "Send to browser":
            #Increment score if solved.
            if responseData["success"] and not synthesis.solverCredited:
                profile = request.user.profile
                profile.correctSynths += 1
                profile.save()
                synthesis.solverCredited = True
                synthesis.save()
            return HttpResponse(json.dumps(responseData))
        elif mode == "Return dictionary":
            return responseData
        
    except StandardError as e:
        responseData = dict()
        responseData["success"] = False
        responseData["molecules"] = [(1, str(e)+traceback.format_exc())]
        responseData["arrows"] = []
        return HttpResponse(json.dumps(responseData))
        
        

@csrf_exempt   
def deleteMolecule(request):    

    #Iteratively delete molecule with id sent in request
    
    
    try:
    
        molIdToDelete = int(request.POST["moleculeID"])
        
        synthesis = request.user.profile.currentSynthesisProblem
        
        markedAny = True
        molIdsToDelete = [molIdToDelete]    #for molecules
        arrIdsToDelete = []                 #for arrows
        
        
        debuggingString = "You said to delete: "+str(molIdToDelete)+"\n"
        
        #While you have recently marked molecules for deletion:
        while markedAny:
            debuggingString += "One loop. \n"
            markedAny = False
            
            #iterate through arrows, checking for any steps with to-delete reactants but not-to-delete products
            for arrowModel in synthesis.arrows.all():
                #if any are found, mark them for deletion
                debuggingString += "Conditions? "+ str(True in [int(arrowModel.pointFrom.id)==int(molId) for molId in molIdsToDelete]) + ", " + str(not (True in [int(arrowModel.pointTo.id)==int(molId) for molId in molIdsToDelete])) + "\n"
                if (True in [int(arrowModel.pointFrom.id)==int(molId) for molId in molIdsToDelete]) and not (True in [int(arrowModel.pointTo.id)==int(molId) for molId in molIdsToDelete]):
                    markedAny = True
                    arrIdsToDelete += [arrowModel.id]
                    molIdsToDelete += [arrowModel.pointTo.id]
                    debuggingString += "Arrow "+str(arrowModel.id)+" with IDs "+str(arrowModel.pointFrom.id)+", "+str(arrowModel.pointTo.id)+" WAS deleted.\n"
                else:
                    debuggingString += "Arrow "+str(arrowModel.id)+" with IDs "+str(arrowModel.pointFrom.id)+", "+str(arrowModel.pointTo.id)+" not deleted.\n"
            
            
        debuggingString += "No more loop! \n"
        
        #Implement: Also iterate through arrows, checking for any steps with to-delete products that are NOT already in the list to delete
        for arrowModel in synthesis.arrows.all():
            if (arrowModel.pointTo.id in molIdsToDelete) and not (arrowModel.id in arrIdsToDelete):
                arrIdsToDelete += [arrowModel.id]
        
        
        #Delete all arrow IDs you found
        for id1 in arrIdsToDelete:
            debuggingString += "Deleted arr: "+str(id1)+"\n"
            try:
                a = models.ArrowModel.objects.get(id=id1)
                synthesis.arrows.remove(a)
                a.delete()
            except:
                raise StandardError(debuggingString)
        
        #Delete all molecule IDs you found
        for id1 in molIdsToDelete:
            debuggingString += "Deleted mol: "+str(id1)+"\n"
            a = models.MoleculeBoxModel.objects.get(id=id1)
            synthesis.molecules.remove(a)
            a.delete()
            
        e = StandardError(debuggingString)
        #raise e
        
        #Return new rendering of problem
        return getSynthesisData(request)

    except BaseException as e:
        responseData = dict()
        responseData["success"] = False
        responseData["molecules"] = [(1, str(e) + traceback.format_exc())]
        responseData["arrows"] = []
        return HttpResponse(json.dumps(responseData))


    
    
def getSolutionData(request):
    solution = request.user.profile.currentSynthesisProblem.solution
    
    #Iterate over all molecules for a specific synthesis
    try:
        moleculesOutput = [ (moleculeBoxModel.id, moleculeBoxModel.svg) 
                        for moleculeBoxModel in solution.molecules.all()]
    except:
        raise Exception("01")
    
    #Iterate over all arrows for a specific synthesis
    arrowsOutput = [ (arrowModel.pointFrom.id, arrowModel.pointTo.id, arrowModel.reagentsHtml)
                     for arrowModel in solution.arrows.all()]
    
    responseData = dict()
    
    responseData["success"] = "solution"
    responseData["molecules"] = moleculesOutput
    responseData["arrows"] = arrowsOutput

    return HttpResponse(json.dumps(responseData))
    
@csrf_exempt
#data: {'molecule1': ui.draggable.attr("id"), 'molecule2': this.attr("id")}
def addMoleculeToMolecule(request):
    #Create a new moleculebox, combining the two inputs into one box, appending it to synthesis.molecules
    #Make sure to check whether they react with each other -- do this via reactionStep
    #Create two new arrows, with empty reagents, correctly attached
    
    if request.method == 'POST':
        try:
            moleculeboxmodel1 = models.MoleculeBoxModel.objects.get(id=request.POST["molecule1"])
            moleculebox1 = moleculeboxmodel1.moleculeBox
            
            moleculeboxmodel2 = models.MoleculeBoxModel.objects.get(id=request.POST["molecule2"])
            moleculebox2 = moleculeboxmodel2.moleculeBox
            
            testStep = ReactionStep(moleculebox1)
            testStep.addMolecule(moleculebox2)
            
            synthesis = request.user.profile.currentSynthesisProblem
            (isTarget, productBox) = testStep.checkStep(synthesis.target.moleculeBox)
            
            productBox.molecules = reactions.removeDuplicates(productBox.molecules)
            
            moleculeboxmodel3 = models.MoleculeBoxModel.create(productBox)
            moleculeboxmodel3.equalsTarget = isTarget
            moleculeboxmodel3.save()
            
            #newPointFrom, newPointTo, newReagentsHtml
            arrow1 = models.ArrowModel.create(moleculeboxmodel1, moleculeboxmodel3, "")
            arrow1.save()
            arrow2 = models.ArrowModel.create(moleculeboxmodel2, moleculeboxmodel3, "")
            arrow2.save()
            
            #Add the arrows and the new moleculebox to the synthesis
            try:
                synthesis.molecules.add(moleculeboxmodel3)
            except:
                raise Exception("02")
            synthesis.arrows.add(arrow1)
            synthesis.arrows.add(arrow2)
            
            
        #For saner debugging
        except StandardError as e:
            responseData = dict()
            responseData["success"] = False
            responseData["arrows"] = []
            responseData["molecules"] = [(1, str(e))]
            return HttpResponse(json.dumps(responseData))
    
    return getSynthesisData(request)

    
    
@csrf_exempt
#data: {'reagents': ui.draggable.attr("reagentString"), 'moleculeOn': this.attr("id")}
def addReagentToMolecule(request):
    #Do this via creating a reactionstep
    #Create a new moleculebox, representing the reaction's products
    #Create a new arrow, containing html of reagents, from reactants to products

    if request.method == 'POST':
        try:
            moleculeboxmodel1 = models.MoleculeBoxModel.objects.get(id=request.POST["moleculeOn"])
            moleculebox1 = moleculeboxmodel1.moleculeBox
            
            reagentString = request.POST["reagents"]
            
            testStep = ReactionStep(moleculebox1)
            testStep.addReagent(parseReagentsString(reagentString))
            
            synthesis = request.user.profile.currentSynthesisProblem
            (isTarget, productBox) = testStep.checkStep(synthesis.target.moleculeBox)
            
            moleculeboxmodel2 = models.MoleculeBoxModel.create(productBox)
            moleculeboxmodel2.equalsTarget = isTarget
            moleculeboxmodel2.save()
        
            arrow = models.ArrowModel.create(moleculeboxmodel1, moleculeboxmodel2, testStep.stringList()[:-2])
            arrow.save()
            
            try:
                synthesis.molecules.add(moleculeboxmodel2)
            except:
                raise Exception("03")
            synthesis.arrows.add(arrow)
        
        except StandardError as e:
            responseData = dict()
            responseData["success"] = False
            responseData["arrows"] = []
            responseData["molecules"] = [(1, str(e))]
            return HttpResponse(json.dumps(responseData))
    
    return getSynthesisData(request)

@csrf_exempt
def askForHelp(request):
  try:
    #Gets called when a user asks for help.
    #If there are no HelpWaitingList's, make one.
    if len(models.HelpWaitingList.objects.all())==0:
        waitingList = models.HelpWaitingList.create()
        waitingList.save()
    else:
        waitingList = models.HelpWaitingList.objects.all()[0]
    #Add this user to the waiting list.
    if len(models.WaitTimer.objects.filter(user=request.user))==0:
        memb = models.WaitTimer(helpWaitingList=waitingList, user=request.user)
        memb.save()
    return HttpResponse(len(waitingList.users.all()))
  except StandardError as e:
    return HttpResponse(str(e))

@csrf_exempt
def helpeeWaitPoll(request):
  try:
    #Gets called every couple of seconds by clients waiting for help.
    out = dict()
    #Look for a chat containing this user.
    chats = models.ChatPair.objects.filter(helpee=request.user)
    if len(chats) == 0 or\
    sum([(timezone.now()-chat.initTime).total_seconds()<TIMEOUT for chat in chats]) == 0:
        #Nope, keep waiting.
        out['success'] = False
        if len(request.user.helpwaitinglist_set.all()) == 0:
            #Oops, the user somehow got wiped from the waiting list.  Better re-add him.
            waitingList = models.HelpWaitingList.objects.all()[0]
            memb = models.WaitTimer(helpWaitingList=waitingList, user=request.user)
            memb.save()
        out['queueSize'] = len(request.user.helpwaitinglist_set.all()[0].users.all())
        #Save this user's waitTimer, so that the lastCheck time updates.
        request.user.waittimer_set.all()[0].save()
        return HttpResponse(json.dumps(out))
    else:
        out['success'] = True
        #Look through all the ChatPairs we found, and delete all the old ones.
        chats = list(chats)
        chats.sort(key=lambda x: x.initTime, reverse=True)
        thisChat = chats[0]
        for i in xrange(1, len(chats)):
            models.ChatPair.delete(chats[i])
        out['helper'] = thisChat.helper.username
        out['chatPK'] = thisChat.pk
        return HttpResponse(json.dumps(out))
  except StandardError as e:
    return HttpResponse(str(e))
    
@csrf_exempt
def volunteerToHelp(request, debug=""):
    #Renders an "instructor lobby" where users who want to help wait for users
    #to need help.
    #Static.
    return render(request, "volunteer.html", {'debug': debug})


@csrf_exempt    
def helperWaitPoll(request):
  try:
    #If there was a POST, that means the helper clicked on a helpee's username.
    #Initiate chat.
    if request.method == 'POST':
        #Get the helpee's username
        helpee = User.objects.get(username=request.POST['username'])
        #Remove the helpee from the queue.
        try:
            models.WaitTimer.delete(helpee.waittimer_set.all()[0])
        except:
            #If this doesn't work, either the helpee has stopped asking for help,
            #or some other helper picked the helpee first.  Either way, go back to
            #the helper lobby and display an error.
            return volunteerToHelp(request, debug="Sorry, that user no longer wants help.")
        #Make a new ChatPair.
        chat = models.ChatPair.create(helpee=helpee, helper=request.user)
        chat.save()
        #Get the helpee's synthesis data.
        out = getSynthesisData(None, synthesis=helpee.profile.currentSynthesisProblem)
        out['target'] = helpee.profile.currentSynthesisProblem.target.svg #Need to send target separately.
        out['helpee'] = helpee.username
        out['chatPK'] = chat.pk
        return HttpResponse(json.dumps(out))
    #Otherwise, we are just updating the list of people who need help.
    #Generates a list of people who need help.
    thisQueue = models.HelpWaitingList.objects.all()[0]
    waitingUsers = ""
    for x in thisQueue.users.all():
        #Did the user poll recently?  If not, he probably closed his browser window.
        #We should remove him from the queue.
        if (timezone.now() - x.waittimer_set.all()[0].lastCheck).total_seconds() > TIMEOUT:
            #Too long.  Remove.
            models.WaitTimer.delete(x.waittimer_set.all()[0])
        else:
            waitingUsers += x.username + " "
    return HttpResponse(waitingUsers)
  except StandardError as e:
    return HttpResponse(str(e))

@csrf_exempt    
def helpeeChatPoll(request):
  try:
    #Is called every couple of seconds by the helpee's browser.
    #Looks for new ChatLines, also updates its helpeeLastCheck.
    if request.method!='POST':
        #This should never happen.
        return
    out = dict()
    pk = int(request.POST['PK'])
    chat = models.ChatPair.objects.get(pk=pk)
    if 'message' in request.POST:
        #A new message was sent.  Add it to the ChatPair.
        msg = escape(request.POST['message'])
        newLine = models.ChatLine.create(originator=request.user, content=msg)
        newLine.save()
        chat.chatRecord.add(newLine)
    
    #Check for helper's messages (regardless of whether a message
    #was sent).
    #Update the last check
    chat.helpeeLastCheck = timezone.now()
    chat.save()
    if (timezone.now() - chat.helperLastCheck).total_seconds() > TIMEOUT:
        #The helper is missing.  Terminate chat.
        out['open'] = False
        models.ChatPair.delete(chat)
        return HttpResponse(json.dumps(out))
    out['open'] = True
    newLines = list(chat.chatRecord.filter(helpeeSeen = False))
    newLines.sort(key=lambda x: x.postTime)
    out['length'] = len(newLines)
    for i in xrange(len(newLines)):
        out[i] = newLines[i].originator.username+ ": " +newLines[i].content
        newLines[i].helpeeSeen = True
        newLines[i].save()
    return HttpResponse(json.dumps(out))
  except BaseException as e:
    return HttpResponse(str(e))
    
@csrf_exempt    
def helperChatPoll(request):
  try:
    #Is called every couple of seconds by the helper's browser.
    #Looks for new ChatLines, also updates its helperLastCheck.
    if request.method!='POST':
        #This should never happen.
        return

    pk = int(request.POST['PK'])

    chat = models.ChatPair.objects.get(pk=pk)
    synthesis = chat.helpee.profile.currentSynthesisProblem
    out = getSynthesisData(None, synthesis)
    if 'message' in request.POST:
        #A new message was sent.  Add it to the ChatPair.
        msg = escape(request.POST['message'])
        newLine = models.ChatLine.create(originator=request.user, content=msg)
        newLine.save()
        chat.chatRecord.add(newLine)

    #Check for helper's messages (regardless of whether a message
    #was sent).
    #Update the last check
    chat.helperLastCheck = timezone.now()
    chat.save()
    if (timezone.now() - chat.helpeeLastCheck).total_seconds() > TIMEOUT:
        #The helper is missing.  Terminate chat.
        out['open'] = False
        models.ChatPair.delete(chat)
        return HttpResponse(json.dumps(out))
    out['open'] = True
    #Send all the new chat lines.
    newLines = list(chat.chatRecord.filter(helperSeen = False))
    newLines.sort(key=lambda x: x.postTime)
    out['length'] = len(newLines)
    for i in xrange(len(newLines)):
        out[i] = newLines[i].originator.username+ ": " + newLines[i].content
        newLines[i].helperSeen = True
        newLines[i].save()
    #If the synthesis was just solved, give the helper credit
    if out["success"] and not synthesis.helperCredited:
        profile = request.user.profile
        profile.assists += 1
        profile.save()
        synthesis.helperCredited = True
        synthesis.save()
    return HttpResponse(json.dumps(out))
  except  BaseException as e:
    return HttpResponse(str(e))
    
@csrf_exempt 
def saveProblem(request):
    #Saves the user's current problem for eternity (sort of).  Sends back an id that
    #can be used to access the problem.
    profile = request.user.profile
    profile.currentSynthesisProblem.retain = True
    profile.currentSynthesisProblem.save()
    return HttpResponse(profile.currentSynthesisProblem.pk)
    
    
    
    
@csrf_exempt 
def renderProblem(request):
    if 'synthesis_resume' in request.POST:
        return renderOldSynthesis(request)
    if 'synthesis_new' in request.POST:
        return renderSynthesis(request)
    if 'synthesis_tutorial' in request.POST:
        return renderSynthesis(request, tutorial=True)
    if 'namereagent_resume' in request.POST:
        return renderOldNameReagent(request)
    if 'namereagent_new' in request.POST:
        return renderNameReagent(request)
    if 'namereagent_tutorial' in request.POST:
        return renderNameReagent(request, tutorial=True)
    else:
        return renderSynthesis(request)
    
    
    
@csrf_exempt
def renderFaq(request):
    return render(request, 'faq.html', {'name': request.user.username})

@csrf_exempt
def renderReactions(request):
    return render(request, 'reactions.html', {'name': request.user.username})
