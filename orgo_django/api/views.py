from django.shortcuts import get_object_or_404, render
from django.http import HttpResponseRedirect, HttpResponse
from django.core.urlresolvers import reverse
from engine.renderSVG import render as svg_render
from engine.randomGenerator import randomStart
from engine.toMolecule import moleculify
from engine.toSmiles import smilesify, to_canonical
import api.engine.reactions

from api.models import Property, Reagent, Reaction, ReagentSet
import json
# Create your views here.

## url(r'^$', views.index, name='index'),
def index(request):
    context = {}
    return render(request, 'api/index.html', context)

## This is so that I can test that toMolecule is working right.
## It's for testing only, and therefore temporary until
## toMolecule has been debugged.
def test_smiles_to_molecule_and_back(request, smiles):
    hydrogens = request.GET.get('hydrogens', False)
    smiles = smiles.replace('/#','#')
    smiles = smiles.replace('~','#')
    molecule = moleculify(smiles)
    smiles2 = '.'.join(smilesify(molecule))
    return HttpResponse('<br />'.join([ smiles, svg_render(smiles, hydrogens), 
                                        smiles2, svg_render(smiles2, hydrogens) ] ))

#/api/findReactions?reagents=[44,2]

#/api/findReagents?text="HBr"

#/api/react?reaction=123&molecules=["SM1","SM2"]

def get_reaction_data(reaction):
    '''Helper function that extracts data from a Reaction query object'''
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
    '''Helper function that extracts data from a Reagent query object'''
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
    '''
    query: a string with the input
    Helper function. Parses the string into a python list of strings
    '''
    if not isinstance(query, (str,unicode) ):
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
    #     #input elements need quotes to convert into strings (eg [aprotic, halide])
    #     #TODO: implement a way to convert the example into ["aprotic","halide"]
    #     msg = "Make sure to use double quotes around string items in your list"
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
    '''
    Creates a JSON list with data objects for each and every reaction
    '''
    reaction_list = Reaction.objects.all().select_related('name', 'id', 'process_function').prefetch_related('reagent_set')
    reaction_data = []
    for reaction in reaction_list:
        reaction_data.append(get_reaction_data(reaction))
    return HttpResponse(json.dumps(reaction_data))

#List basic info about a single reaction, id#123 in the database [as JSON]
def get_reaction(request, id):
    '''Creates a JSON object with data for the reaction specified by id'''
    reaction = Reaction.objects.get(id=id)
    reaction_data = get_reaction_data(reaction)
    return HttpResponse(json.dumps(reaction_data))

#If the user entered in [list of reagents, by id], what reaction(s) do I get?
def find_reactions(request):
    if request.method == "GET":
        reagent_id_list = parse_input(request.GET.get('reagents', ''))
        reagent_set_list = ReagentSet.objects.all()
        #Recursively filter for each reagent to get only reagent sets that have ALL the reagents
        for reagent_id in reagent_id_list:
            reagent_set_list = reagent_set_list.filter(reagents__id=int(reagent_id))

        reaction_list = reagent_set_list.reactions.all()

        reaction_data = []
        for reaction in reaction_list:
            reaction_data.append(get_reaction_data(reaction))
    return HttpResponse(json.dumps(reaction_data))

#List links to all reagents
def all_reagents(request):
    reagent_list = Reagent.objects.all()
    reagent_data = []
    for reagent in reagent_list:
      reagent_data.append(get_reagent_data(reagent))
    return HttpResponse(json.dumps(reagent_data))

#List basic info about a single reagent, id#123 in the database [as JSON]
def get_reagent(request, id):
    reagent = Reagent.objects.get(id=id)
    reagent_data = get_reagent_data(reagent)
    return HttpResponse(json.dumps(reagent_data))

#use to get a valid reagent by name. throws error otherwise
def get_valid_reagent(request):
    name = request.GET.get("name", '')
    name = name.strip('"').strip("'")
    if len(name) == 0:
        raise error("No name provided")
    reagent = Reagent.objects.get(name__iexact=name)
    return HttpResponse(reagent.name)

#If the user entered in name=input or properties=input, what reagent(s) do I get?
def find_reagents(request):
    if request.method == "GET":
        reagent_name = request.GET.get('name', '')
        reagent_name = reagent_name.strip('"').strip("'")

        reagent_props = parse_input(request.GET.get('properties',''))

        reagent_list = Reagent.objects.all()

        if reagent_name is not '':
            reagent_list = reagent_list.filter(name__istartswith=reagent_name) #TODO: Change once name is changed to StringListField
        
        #Recursively filter for each property to get only reagents that have ALL the properties
        for prop_name in reagent_props:
            #iexact for case-insensitive matches (eg Aprotic = aprotic = apROtiC)
            reagent_list = reagent_list.filter(properties__name__iexact=prop_name)

        reagent_list = reagent_list.select_related('name', 'id', 'diagram_name', 'smiles', 'is_solvent').prefetch_related('properties')
        reagent_data = []
        for reagent in reagent_list:
            reagent_data.append(get_reagent_data(reagent))
    return HttpResponse(json.dumps(reagent_data))


#what if the SMILES are different but represent the same molecule?
def check_if_equal(request):
    '''Return true if two SMILES represent the same molecule'''
    mol1 = request.GET.get('mol1', None)
    mol2 = request.GET.get('mol2', None)
    mol1 = to_canonical(mol1)
    mol2 = to_canonical(mol2)
    return HttpResponse(mol1 == mol2)

#React molecule(s) (by SMILES) with a particular reaction, and return the result (as SMILES)
def react(request):
    #pseudocode/almost code but it wouldn't actually work
    # TODO: Make this work!
    reactionID = request.GET.get('reaction', None)
    reactants = moleculify(request.GET.get('reactants', None))
    function_name = Reaction.objects.get(id = reactionID).process_function
    #raise StandardError(function_name)
    reaction = getattr(api.engine.reactions, function_name)
    products = reaction(reactants)
    return HttpResponse(smilesify(products))

#Render a molecule (convert SMILES to SVG)
def render_SVG(request):
    smiles = request.GET.get('molecule', None) # change to error raise later
    hydrogens = request.GET.get('hydrogens', False)
    return HttpResponse(svg_render(smiles, hydrogens)) # hydrogens default to false

#Return a randomly-generated molecule (output a SMILES)
def random_gen_smiles(request):
    mol, _, _ = randomStart()
    return HttpResponse(smilesify(mol))

#Render a randomly-generated molecule (output a SVG)
def random_gen_SVG(request):
    mol, _, _ = randomStart()
    return HttpResponse(svg_render(smilesify(mol)))
