from django.shortcuts import render
from django.shortcuts import get_object_or_404, render
from django.http import HttpResponseRedirect, HttpResponse
from django.core.urlresolvers import reverse

from random import randrange

from api.models import Reagent
from app.models import SingleStepProblem
NUM_OPTIONS = 4

# Create your views here.

## url(r'^$', views.index, name='index'),
def index(request):
    context = {}
    return render(request, 'app/index.html', context)

def synthesis(request):
    return HttpResponse("hi")
    # pseudocode:
    # while predict_products(request)!=desiredProduct: (figure out how to desiredProduct)
    #    add_request
    #    get_new_request

def single_step(request, id):
    context = {}
    problem = SingleStepProblem.objects.get(id=id)
    context["reactantSmiles"] = problem.reactantSmiles
    context["productSmiles"] = problem.productSmiles
    context["correctAnswer"] = problem.correctAnswer
    
    context["NUM_OPTIONS"] = NUM_OPTIONS
    numReagents = len(Reagent.objects.all())
    options = []
    for i in xrange(NUM_OPTIONS - 1):
        reagent = Reagent.objects.get(id=randrange(numReagents) + 1) # can break
        options.append(reagent)
    context["incorrectAnswers"] = options
    context["incorrectAnswer"] = "wrong answer"
    return render(request, 'app/singleStep.html', context)

def predict_products(request):
    ## NOT PERFECTLY DOABLE YET
    ## possibleReactions = findReactions(request)
    ## products = react(request)
    return True #### TODO: figure out what in the world is going on
