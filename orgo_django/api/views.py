"""
api/views

Contains API views.
"""

from django.shortcuts import render
from django.http import HttpResponse, Http404
from engine.renderSVG import render as svg_render
from engine.randomGenerator import random_molecule
from engine.toMolecule import moleculify
from engine.toSmiles import smilesify, to_canonical
import api.engine.reaction_functions

from api.models import Reagent, Reaction, ReagentSet
import json
# Create your views here.

##### HELPER FUNCTIONS #####
def JsonResponse(data):
    return HttpResponse(json.dumps(data, sort_keys=True, indent=4),
                        content_type="application/json")


##### VIEWS #####

## url(r'^$', views.index, name='index'),
def index(request):
    "Default index page"
    context = {}
    return render(request, 'api/index.html', context)

## This is so that I can test that toMolecule is working right.
## It's for testing only, and therefore temporary until
## toMolecule has been debugged.
def test_parsers(request, smiles):
    "View for visually testing toMolecule/toSmiles"
    hydrogens = request.GET.get('hydrogens', False)
    smiles = smiles.replace('/#', '#')
    smiles = smiles.replace('~', '#')
    molecule = moleculify(smiles)
    smiles2 = smilesify(molecule, canonical=False)
    # print "Smiles produced: { %s }" % smiles2
    return HttpResponse('<br />'.join([
        "<h2>Smiles-to-Molecule Testing Viewer</h2>",
        "<b>Original:</b> %s" % smiles,
        "<b>Canonical:</b> %s" % to_canonical(smiles),
        svg_render(to_canonical(smiles), hydrogens),
        "<b>After two passes:</b> %s" % smiles2,
        "<b>Canonical:</b> %s" % to_canonical(smiles2),
        svg_render(to_canonical(smiles2), hydrogens),
        "<h1>%s</h1>" % str((to_canonical(smiles) == to_canonical(smiles2))),
    ]))

#/api/findReactions?reagents=[44,2]

#/api/findReagents?text="HBr"

#/api/react?reaction=123&molecules=["SM1","SM2"]

def get_reaction_data(reaction):
    """Helper function that extracts data from a Reaction query object"""
    solvent = reaction.reagent_set.solvent
    if solvent is not None:
        solvent = solvent.id
    data = {
        "id": reaction.id,
        "name": reaction.name,
        "process_function": reaction.process_function,
        "reagents": [r.id for r in reaction.reagent_set.reagents.all()],
        "solvent": solvent,
        "solvent_properties": [
            prop.name for prop in reaction.reagent_set.solvent_properties.all()
        ],
    }
    return data

def get_reagent_data(reagent):
    """Helper function that extracts data from a Reagent query object"""
    data = {
        "id": reagent.id,
        "name": reagent.name, #TODO: Change if name becomes a StringListField
        "is_solvent": reagent.is_solvent,
        "diagram_name": reagent.diagram_name,
        "smiles": reagent.smiles,
        "properties": [prop.name for prop in reagent.properties.all()]
    }
    return data

#allows the query input to be in several forms: [x,y], ["x,y"], x,y, etc.
def parse_input(query):
    """
    query: a string with the input
    Helper function. Parses the string into a python list of strings
    """
    if not isinstance(query, (str, unicode)):
        print "Check code. Should not be parsing non-strings"
        return query
    if query == '':
        return  []
    temp = query.strip('()[]{}').split(',')
    parsed_query = [i.strip('"').strip("'") for i in temp]
    return parsed_query

    # try:
    #     #convert from str/unicode to an actual python list (or string)
    #     parsed_query = json.loads(query)
    # except ValueError:
    #     #input elements need quotes to convert into str (eg [aprotic, halide])
    #     #TODO: implement a way to convert example into ["aprotic","halide"]
    #     msg = "Make sure to use double quotes around string items in list"
    #     raise ValueError(msg)
    # except TypeError:
    #     msg = "Cannot parse input data (list). Provide input as valid list"
    #     raise TypeError(msg)
    # except Exception as e:
    #     msg = "Unexpected error: " + type(e)
    #     print msg #this print probably doesn't work
    #     raise
    # else:
    #     return parsed_query

#List reaction data for ALL reactions
def all_reactions(request):
    """
    Creates a JSON list with data objects for each and every reaction
    """
    reactions = Reaction.objects.all()\
        .select_related('name', 'id', 'process_function')\
        .prefetch_related('reagent_set')
    reaction_data = []
    for reaction in reactions:
        reaction_data.append(get_reaction_data(reaction))
    return JsonResponse(reaction_data)

#List basic info about a single reaction, id#123 in the database [as JSON]
def get_reaction(request, pk):
    """Creates a JSON object with data for the reaction specified by pk"""
    reaction = Reaction.objects.get(id=pk)
    reaction_data = get_reaction_data(reaction)
    return JsonResponse(reaction_data)

#If the user entered in [list of reagents, by id], what reaction(s) do I get?
def find_reactions(request):
    ##TODO
    return HttpResponse("Currently buggy -- not yet fixed. TODO.")
    if request.method == "GET":
        reagent_ids = parse_input(request.GET.get('reagents', ''))
        reagent_sets = ReagentSet.objects.all()
        #Recursively filter to get only reagent sets that have ALL the reagents
        for reagent_id in reagent_ids:
            reagent_sets = reagent_sets.filter(reagents__id=int(reagent_id))

        reactions = reagent_sets.reactions.all()

        reaction_data = []
        for reaction in reactions:
            reaction_data.append(get_reaction_data(reaction))
    return JsonResponse(reaction_data)

#List links to all reagents
def all_reagents(request):
    reagents = Reagent.objects.all()
    reagent_data = []
    for reagent in reagents:
        reagent_data.append(get_reagent_data(reagent))
    return JsonResponse(reagent_data)

#List basic info about a single reagent, id#123 in the database [as JSON]
def get_reagent(request, pk):
    reagent = Reagent.objects.get(id=pk)
    reagent_data = get_reagent_data(reagent)
    return JsonResponse(reagent_data)

#use to get a valid reagent by name. throws error otherwise
# def get_valid_reagent(request):
#     name = request.GET.get("name", '')
#     name = name.strip('"').strip("'")
#     if len(name) == 0:
#         raise StandardError("No name provided")
#     reagent = Reagent.objects.get(name__iexact=name)
#     return HttpResponse(reagent.name)

#If the user enters in name=input or properties=input, what reagent(s) do I get?
def find_reagents(request):
    if request.method == "GET":
        reagent_name = request.GET.get('name', '')
        reagent_name = reagent_name.strip('"').strip("'")
        reagent_props = parse_input(request.GET.get('properties', ''))
        reagents = Reagent.objects.all()

        if reagent_name is not '':
            reagents = reagents.filter(name__istartswith=reagent_name)
                # ^ TODO: Change once name is changed to StringListField
        
        #Recursively filterto get only reagents that have ALL the properties
        for prop_name in reagent_props:
            #iexact for case-insensitive (eg Aprotic = aprotic = apROtiC)
            reagents = reagents.filter(properties__name__iexact=prop_name)

        reagents = reagents.select_related(
            'name', 'id', 'diagram_name',
            'smiles', 'is_solvent'
        ).prefetch_related('properties')
        
        reagent_data = []
        for reagent in reagents:
            reagent_data.append(get_reagent_data(reagent))
    return JsonResponse(reagent_data)


#what if the SMILES are different but represent the same molecule?
def check_if_equal(request):
    """Return true if two SMILES represent the same molecule"""
    mol1 = request.GET.get('mol1', None)
    mol2 = request.GET.get('mol2', None)
    mol1 = to_canonical(mol1)
    mol2 = to_canonical(mol2)
    return HttpResponse(mol1 == mol2)

def react(request):
    """
    React molecule(s) (by SMILES) with a particular reaction, and return
    the result (as SMILES)
    """
    #pseudocode/almost code but it wouldn't actually work
    # TODO: Make this work!
    reactionID = request.GET.get('reaction', None)
    reactants = json.loads(request.GET.get('reactants', None))
    if reactants == None:
        return ""
    reactants = moleculify(reactants)
    try:
        function_name = Reaction.objects.get(id=reactionID).process_function
    except:
        raise Http404
    reaction = getattr(api.engine.reaction_functions, function_name)
    products = reaction(reactants)
    return HttpResponse(smilesify(products))

#Render a molecule (convert SMILES to SVG)
def render_SVG(request):
    smiles = request.GET.get('molecule', None) # change to error raise later TODO
    if smiles == None:
        smiles = request.GET.get('mol', None)
    return HttpResponse(svg_render(smiles))

#Return a randomly-generated molecule (output a SMILES)
def random_gen_smiles(request):
    mol = random_molecule()
    print "Made! Smilesifying..."
    return HttpResponse(smilesify(mol))

#Render a randomly-generated molecule (output a SVG)
def random_gen_SVG(request):
    mol = random_molecule()

    return HttpResponse(svg_render(smilesify(mol)))

def to_canonical_view(request):
    smiles = request.GET.get('molecule', None) # change to error raise later TODO
    if smiles == None:
        smiles = request.GET.get('mol', None)
    return HttpResponse(to_canonical(smiles))