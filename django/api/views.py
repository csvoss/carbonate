"""
api/views

Contains API views.
"""

from django.shortcuts import render
from django.http import HttpResponse, HttpResponseBadRequest, HttpResponseNotFound
from engine.renderSVG import render as svg_render
from engine.randomGenerator import random_molecule
from engine.toMolecule import moleculify
from engine.toSmiles import smilesify, to_canonical
import api.engine.reaction_functions

from api.models import Property, ReagentName, Reagent, Reaction, ReagentSet
import json

############################
##### HELPER FUNCTIONS #####
############################

def JsonResponse(data):
    """
    data :: dict or list.
    return :: HttpResponse. Pretty-printed and JSON-formatted.
    """
    return HttpResponse(json.dumps(data, sort_keys=True, indent=4),
                        content_type="application/json")

def SvgResponse(data):
    """
    data :: SVG string.
    return :: HttpResponse. SVG with proper content-type header.
    """
    return HttpResponse(data, content_type="image/svg")

def RequiredQueryParamsResponse(param):
    return HttpResponseBadRequest('Must specify query parameter "%s" in your request.' % str(param))

def ModelNotFoundResponse(model, pk):
    return HttpResponseNotFound('No %s found by that identifier: %s' % (model, str(pk)))

def get_reaction_data(reaction):
    """
    Helper function that extracts data from a Reaction query object.
    reaction :: models.Reaction.
    return :: dict.
    """
    solvent = reaction.reagent_set.solvent
    if solvent is not None:
        solvent = solvent.id
    data = {
        "id": reaction.id,
        "name": reaction.name,
        "process_function": reaction.process_function,
        "reagents": [r.id for r in reaction.reagent_set.reagents.all()],
        "solvent": solvent,
        "solvent_properties": [prop.name for prop in reaction.reagent_set.solvent_properties.all()],
    }
    return data

def get_reagent_data(reagent):
    """
    Helper function that extracts data from a Reagent query object.
    reagent :: models.Reagent.
    return :: dict.
    """
    data = {
        "id": reagent.id,
        "names": [reagent_name.name for reagent_name in reagent.names.all()],
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
    if not isinstance(query, basestring):
        print "Check code. Should not be parsing non-strings"
        return query
    if query == '':
        return  []
    temp = query.strip('()[]{}').split(',')
    parsed_query = [i.strip('"').strip("'") for i in temp]
    return parsed_query

#################
##### VIEWS #####
#################

## url(r'^$', views.index, name='index'),
def index(request):
    "Default index page"
    return render(request, 'api/index.html')

def test_parsers(request, smiles):
    """
    This is so that I can test that toMolecule is working right.
    It's for testing only, and therefore temporary until
    toMolecule has been debugged.
    """
    smiles = smiles.replace('/#', '#')
    smiles = smiles.replace('~', '#')
    molecule = moleculify(smiles)
    smiles2 = smilesify(molecule, canonical=False)
    return HttpResponse('<br />'.join([
        "<h2>Smiles-to-Molecule Testing Viewer</h2>",
        "<b>Original:</b> %s" % smiles,
        "<b>Canonical:</b> %s" % to_canonical(smiles),
        svg_render(to_canonical(smiles)),
        "<b>After two passes:</b> %s" % smiles2,
        "<b>Canonical:</b> %s" % to_canonical(smiles2),
        svg_render(to_canonical(smiles2)),
        "<h1>%s</h1>" % ("Good!" if to_canonical(smiles) == to_canonical(smiles2) else "Bad..."),
    ]))



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
    try:
        reaction = Reaction.objects.get(id=pk)
    except Reaction.DoesNotExist:
        return ModelNotFoundResponse("reaction", str(pk))
    reaction_data = get_reaction_data(reaction)
    return JsonResponse(reaction_data)

#If the user entered in [list of reagents, by id], what reaction(s) do I get?
def find_reactions(request):
    ##TODO -- seems buggy
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
    try:
        reagent = Reagent.objects.get(id=pk)
    except Reagent.DoesNotExist:
        return ModelNotFoundResponse("reagent", str(pk))
    reagent_data = get_reagent_data(reagent)
    return JsonResponse(reagent_data)

def all_reagent_names(request):
    reagent_data = []
    for rname in ReagentName.objects.all():
        other_names = "Aka: "
        for rn in rname.reagent.names.all():
            if rn.name != rname.name:
                other_names += rn.name + ", "
        other_names = other_names.rstrip(", ")
        name_data = {
            "name": rname.name,
            "id": rname.reagent.id,
            "description": other_names
        }
        reagent_data.append(name_data)
    return JsonResponse(reagent_data)

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


def check_if_equal(request):
    """
    request     mol1 :: str (SMILES).
                mol2 :: str (SMILES).
    Return true if two SMILES represent the same molecule.
    Canonical SMILES is the best!
    """
    mol1 = request.GET.get('mol1', None)
    mol2 = request.GET.get('mol2', None)

    if mol1 == None:
        return RequiredQueryParamsResponse('mol1')
    if mol2 == None:
        return RequiredQueryParamsResponse('mol2')

    mol1_can = to_canonical(mol1)
    mol2_can = to_canonical(mol2)
    if mol1_can == "" or mol2_can == "":
        return JsonResponse(mol1 == mol2) ## Invalid SMILES are not equal
    return JsonResponse(mol1_can == mol2_can)

def react(request):
    """
    React molecule(s) (by SMILES) with a particular reaction, and return
    the result (as SMILES)
    """
    #pseudocode/almost code but it wouldn't actually work
    # TODO: Make this work!
    reactionID = request.GET.get('reaction', None)
    reactant_smi_list = request.GET.get('reactants', None)

    print reactant_smi_list

    if reactionID == None:
        return RequiredQueryParamsResponse('reaction')
    if reactant_smi_list == None:
        return JsonResponse("")

    try:
        reactant_smi_list = json.loads(reactant_smi_list)
    except ValueError:
        return JsonResponse("")

    reactants = moleculify(reactant_smi_list)
    
    try:
        reaction = Reaction.objects.get(id=reactionID)
    except Reaction.DoesNotExist:
        return ModelNotFoundResponse("reaction", str(reactionID))
    
    function_name = reaction.process_function
    reaction_function = getattr(api.engine.reaction_functions, function_name)
    products = reaction_function(reactants)
    return JsonResponse(smilesify(products))

def react_then_SVG(request):
    """
    For testing purposes.
    """
    reactionID = request.GET.get('reaction', None)
    reactant_smi_list = request.GET.get('reactants', None)

    if reactionID == None:
        return RequiredQueryParamsResponse('reaction')
    if reactant_smi_list == None:
        return HttpResponse("")

    try:
        reactant_smi_list = json.loads(reactant_smi_list)
    except ValueError:
        return JsonResponse("")

    reactants = moleculify(reactant_smi_list)
    
    try:
        reaction = Reaction.objects.get(id=reactionID)
    except Reaction.DoesNotExist:
        return ModelNotFoundResponse("reaction", str(reactionID))
    
    function_name = reaction.process_function
    reaction_function = getattr(api.engine.reaction_functions, function_name)
    products = reaction_function(reactants)
    return HttpResponse("<br/".join([
        reaction.name,
        svg_render(to_canonical('.'.join(reactant_smi_list))),
        svg_render(smilesify(products)),
    ]))

def render_SVG(request):
    """
    Render a molecule (convert SMILES to SVG)
    """
    smiles = request.GET.get('molecule', request.GET.get('mol', None))
    if smiles == None:
        return RequiredQueryParamsResponse('mol')
    return SvgResponse(svg_render(smiles))


def random_smiles(request):
    """
    Return a randomly-generated molecule (output a SMILES)
    """
    seed = request.GET.get('seed', None)
    mol = random_molecule(seed=seed)
    smi = smilesify(mol, canonical=True)
    return JsonResponse(smi)

def random_SVG(request):
    """
    Render a randomly-generated molecule (output a SVG)
    """
    seed = request.GET.get('seed', None)
    mol = random_molecule(seed=seed)
    return SvgResponse(svg_render(smilesify(mol)))

def to_canonical_view(request):
    """
    Convert a SMILES string (in the query parameters) to canonical SMILES
    (returned as a JSON-encoded string).
    """
    smiles = request.GET.get('molecule', request.GET.get('mol', None))
    if smiles == None:
        return RequiredQueryParamsResponse('mol')
    return JsonResponse(to_canonical(smiles))

def is_correct_reagent_set(request):
    submitted_reagent_names = request.GET.get('submitted_reagents', [])
    submitted_solvent_id = request.GET.get('solvent_id', None)
    solution_set_ids = request.GET.get('solution_ids', [])

    submitted_reagent_ids = sorted( [ReagentName.objects.get(name=rname).reagent.id for rname in submitted_reagent_names] )

    correct = False
    for set_id in solution_set_ids:
        valid_set = ReagentSet.objects.get(id=set_id)
        reagent_ids = sorted( [r.id for r in valid_set.reagents.all()] )

        solvent_id = valid_set.solvent.id
        if reagent_ids == submitted_reagent_ids and submitted_solvent_id == solvent_id:
            correct = True

    return HttpResponse(correct)

def update_reagent_names(request):
    reagents = Reagent.objects.all()
    createdNames = []
    for rg in reagents:
        name = unicode(rg)
        (obj, created) = ReagentName.objects.get_or_create(name=name, reagent=rg)
        if created:
            createdNames.append(name)
    return JsonResponse(createdNames)
