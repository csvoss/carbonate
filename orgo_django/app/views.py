from django.shortcuts import get_object_or_404, render
from django.http import HttpResponseRedirect, HttpResponse
from django.core.urlresolvers import reverse
import json
from random import randrange
from api.models import Reagent
from api.models import Reaction
from app.models import Synthesis
from app.models import SingleStepProblem
from app.models import PredictProductsProblem

NUM_OPTIONS = 4

# Create your views here.
## url(r'^$', views.index, name='index'),
def index(request):
    context = {}
    return render(request, 'app/index.html', context)

def synthesis(request):
    context = {}
    problem = SingleStepProblem.objects.get(id=id)
    context["reactantSmiles"] = problem.reactantSmiles
    context["productSmiles"] = problem.productSmiles
    context["correctAnswer"] = problem.correctAnswer
    return render(request, 'app/synthesis.html', context)
    # for future reference
    # pseudocode:
    # while predict_products(request)!=desiredProduct: (figure out how to desiredProduct)
    #    add_request
    #    get_new_request
    # return render(request, 'app/synthesis.html', context)

def single_step(request, id):
    context = {}
    problem = SingleStepProblem.objects.get(id=id)
    context["reactantSmiles"] = problem.reactantSmiles
    context["productSmiles"] = problem.productSmiles
    context["correctAnswer"] = problem.correctAnswer
    context["NUM_OPTIONS"] = NUM_OPTIONS
    numReagents = len(Reagent.objects.all())
    options = []
    correctIndex = randrange(NUM_OPTIONS)
    for i in xrange(NUM_OPTIONS):
        if  (i == correctIndex):
            options.append(problem.correctAnswer)
        else:
            reagent = Reagent.objects.get(id=randrange(numReagents) + 1)
            reagentName = reagent.name
            options.append(reagentName)
    context["answers"] = json.dumps(options)
    return render(request, 'app/singleStep.html', context)

def predict_products(request, id):
    context = {}
    problem = PredictProductsProblem.objects.get(id=id)
    context["reactantSmiles"] = problem.reactantSmiles
    context["reagents"] = problem.reagents
    context["correctAnswer"] = problem.correctAnswer
    context["NUM_OPTIONS"] = NUM_OPTIONS
    numReactions = len(Reaction.objects.all())
    options = []
    correctIndex = randrange(NUM_OPTIONS)
    firstOptionAvailable = True
    secondOptionAvailable = True
    for i in xrange(NUM_OPTIONS):
        if  (i == correctIndex):
            options.append(problem.correctAnswer)
        elif (firstOptionAvailable):
            options.append(problem.incorrectAnswer1)
            firstOptionAvailable = False
        elif (secondOptionAvailable):
            options.append(problem.incorrectAnswer2)
            secondOptionAvailable = False
        else:
            options.append(problem.incorrectAnswer3)
            
    context["answers"] = json.dumps(options)
    return render(request, 'app/predictProducts.html', context)
    
    ## MAYBE THIS ACTUALLY WORKS
    ## NOT PERFECTLY DOABLE YET
    ## possibleReactions = findReactions(request)
    ## products = react(request)
    #### TODO: figure out what in the world is going on
