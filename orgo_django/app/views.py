"""
views.py
Convert Django's HTTP requests, routed here by urls.py,
into responses to return to the user.
"""

from django.shortcuts import get_object_or_404, render
from django.http import HttpResponseRedirect, HttpResponse
from django.core.urlresolvers import reverse
import json
from random import randrange, shuffle
from api.models import Property, Reagent, ReagentSet, Reaction

from app.models import Synthesis, SingleStepProblem, SingleStepHardProblem, PredictProductsProblem

NUM_OPTIONS = 4

# Create your views here.
## url(r'^$', views.index, name='index'),
def index(request):
    """
    Return a page for /app.
    request :: HttpRequest
    return :: HttpResponse
    """
    context = {}
    return render(request, 'app/index.html', context)

def synthesis(request, id):
    """
    Create a UI for SynthesisProblem #id.
    request :: HttpRequest
    id :: int
    return :: HttpResponse
    """
    problem = Synthesis.objects.get(id=id)
    context = {
        "reactant_smiles": problem.reactant_smiles,
        "product_smiles": problem.product_smiles,
        "correct_answer": problem.correct_answer,
    }
    return render(request, 'app/synthesis.html', context)
    # for future reference
    # pseudocode:
    # while predict_products(request)!=desiredProduct: (figure out how to desiredProduct)
    #    add_request
    #    get_new_request
    # return render(request, 'app/synthesis.html', context)

def single_step(request, id):
    """
    Create a UI for SingleStepProblem #id.
    request :: HttpRequest
    id :: int
    return :: HttpResponse
    """
    num_reagents = len(Reagent.objects.all())
    correct_index = randrange(NUM_OPTIONS)

    # This way, no reagents will be repeated
    reagents = range(1,num_reagents+1)
    shuffle(reagents)
    reagent_choices = reagents[:NUM_OPTIONS-1] 

    problem = SingleStepProblem.objects.get(id=id)
    
    options = [Reagent.objects.get(id=i).name for i in reagent_choices]
    options.insert(correct_index, problem.correct_answer)

    context = {
        "reactant_smiles": problem.reactant_smiles,
        "product_smiles": problem.product_smiles,
        "correct_answer": problem.correct_answer,
        "NUM_OPTIONS": NUM_OPTIONS,
        "answers": json.dumps(options)
    }
    return render(request, 'app/singleStep.html', context)

def single_step_hard(request, id):
    """
    Create a UI for SingleStepHardProblem #id.
    request :: HttpRequest
    id :: int
    return :: HttpResponse
    """
    problem = SingleStepHardProblem.objects.get(id=id)
    # TODO: account for solvent being reagent or properties
    # solvent = problem.answer.solvent.name
    context = {
        'reactant': problem.reactant_smiles,
        'product' : problem.product_smiles,
        'reagents': json.dumps([reagent.name for reagent in problem.answer.reagents.all()]),
        # 'solvent' : solvent,
    }
    return render(request, 'app/SingleStepHard.html', context)

def predict_products(request, id):
    """
    Create a UI for PredictProductsProblem #id.
    request :: HttpRequest
    id :: int
    return :: HttpResponse
    """
    problem = PredictProductsProblem.objects.get(id=id)
    correct_index = randrange(NUM_OPTIONS)
    options = [
        problem.incorrect_answer1,
        problem.incorrect_answer2,
        problem.incorrect_answer3,
    ]
    shuffle(options)
    options.insert(correct_index, problem.correct_answer)
    context = {
        "reactant_smiles": problem.reactant_smiles,
        "reagents": problem.reagents,
        "correct_answer": problem.correct_answer ,
        "NUM_OPTIONS": NUM_OPTIONS,
        "answers": json.dumps(options),
    }
            
    return render(request, 'app/predictProducts.html', context)

    ## MAYBE THIS ACTUALLY WORKS
    ## NOT PERFECTLY DOABLE YET
    ## possibleReactions = findReactions(request)
    ## products = react(request)
    #### TODO: figure out what in the world is going on
