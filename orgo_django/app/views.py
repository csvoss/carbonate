from django.shortcuts import get_object_or_404, render
from django.http import HttpResponseRedirect, HttpResponse
from django.core.urlresolvers import reverse
import json
from random import randrange
from api.models import Property, Reagent, ReagentSet, Reaction

from app.models import Synthesis, SingleStepProblem, SingleStepHardProblem, PredictProductsProblem

NUM_OPTIONS = 4

# Create your views here.
## url(r'^$', views.index, name='index'),
def index(request):
    context = {}
    return render(request, 'app/index.html', context)

def synthesis(request, id):
    context = {}
    problem = SingleStepProblem.objects.get(id=id)
    context["reactant_smiles"] = problem.reactant_smiles
    context["product_smiles"] = problem.product_smiles
    context["correct_answer"] = problem.correct_answer 
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
    context["reactant_smiles"] = problem.reactant_smiles
    context["product_smiles"] = problem.product_smiles
    context["correct_answer"] = problem.correct_answer 
    context["NUM_OPTIONS"] = NUM_OPTIONS
    numReagents = len(Reagent.objects.all())
    options = []
    correctIndex = randrange(NUM_OPTIONS)
    for i in xrange(NUM_OPTIONS):
        if  (i == correctIndex):
            options.append(problem.correct_answer)
        else:
            reagent = Reagent.objects.get(id=randrange(numReagents) + 1)
            reagentName = reagent.name
            options.append(reagentName)
    context["answers"] = json.dumps(options)
    return render(request, 'app/singleStep.html', context)

def single_step_hard(request, id):
    problem = SingleStepHardProblem.objects.get(id=id)
    # solvent = problem.answer.solvent.name
    return render(request, 'app/SingleStepHard.html', {
        'reactant': problem.reactant_smiles,
        'product' : problem.product_smiles,
        'reagents': [reagent.name for reagent in problem.answer.reagents.all()],
        # 'solvent' : solvent,
        })

def predict_products(request, id):
    context = {}
    problem = PredictProductsProblem.objects.get(id=id)
    context["reactant_smiles"] = problem.reactant_smiles
    context["reagents"] = problem.reagents
    context["correct_answer"] = problem.correct_answer 
    context["NUM_OPTIONS"] = NUM_OPTIONS
    numReactions = len(Reaction.objects.all())
    options = []
    correctIndex = randrange(NUM_OPTIONS)
    firstOptionAvailable = True
    secondOptionAvailable = True
    for i in xrange(NUM_OPTIONS):
        if  (i == correctIndex):
            options.append(problem.correct_answer  )
        elif (firstOptionAvailable):
            options.append(problem.incorrect_answer1)
            firstOptionAvailable = False
        elif (secondOptionAvailable):
            options.append(problem.incorrect_answer2)
            secondOptionAvailable = False
        else:
            options.append(problem.incorrect_answer3)
            
    context["answers"] = json.dumps(options)
    return render(request, 'app/predictProducts.html', context)
    
    ## MAYBE THIS ACTUALLY WORKS
    ## NOT PERFECTLY DOABLE YET
    ## possibleReactions = findReactions(request)
    ## products = react(request)
    #### TODO: figure out what in the world is going on
