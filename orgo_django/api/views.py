from django.shortcuts import get_object_or_404, render
from django.http import HttpResponseRedirect, HttpResponse
from django.core.urlresolvers import reverse
from engine.renderSVG import render as svg_render
from engine.randomGenerator import randomStart
from engine.toMolecule import moleculify
from engine.toSmiles import smilesify, to_canonical
import api.engine.reactions

from api.models import Property, Reagent, Reaction
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
    smiles2 = smilesify(molecule)
    return HttpResponse('<br />'.join([ smiles, svg_render(smiles, hydrogens), 
                                        smiles2, svg_render(smiles2, hydrogens) ] ))

#/api/findReactions?reagents=[44,2]

#/api/findReagents?text="HBr"

#/api/react?reaction=123&molecules=["SM1","SM2"]


def get_reaction_data(reaction):
    '''Helper function that extracts data from a Reaction query object'''
    data = {
        "id": reaction.id,
        "name": reaction.name,
        "process_function": reaction.process_function,
        "reagents": [r.id for r in reaction.reagents.all()],
        "solvent": reaction.solvent.id,
        "solvent_properties": [prop.name for prop in reaction.solvent_properties.all()],
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

#List reaction data for ALL reactions
def reactions(request):
    '''
    Creates a JSON list with data objects for each and every reaction
    '''
    reaction_list = Reaction.objects.all().select_related('name', 'id', 'process_function', 'solvent').prefetch_related('reagents', 'solvent_properties')
    reaction_data = []
    for reaction in reaction_list:
        reaction_data.append(get_reaction_data(reaction))
    return HttpResponse(json.dumps(reaction_data))

#List basic info about a single reaction, id#123 in the database [as JSON]
def reaction(request, id):
    '''Creates a JSON object with data for the reaction specified by id'''
    reaction = Reaction.objects.get(id=id)
    reaction_data = get_reaction_data(reaction)
    return HttpResponse(json.dumps(reaction_data))

#If the user entered in [list of reagents, by id], what reaction(s) do I get?
def find_reactions(request):
    if request.method == "GET":
        reagent_id_list = request.GET.get('reagents', None)
        if isinstance(reagent_id_list, (str, unicode)):
            #convert from str/unicode to an actual python list
            reagent_id_list = json.loads(reagent_id_list)
        elif isinstance(reagent_id_list, list):
            #keep the list as is if it's already a list for some reason
            pass
        elif reagent_id_list == None:
            #input is empty, show all reactions
            return reactions(request)
        else:
            #invalid input?
            msg = "Invalid Input. Provide reagents as a list of ids"
            return HttpResponse(msg)
        reaction_list = Reaction.objects.all()
        #Recursively filter for each reagent to get only reactions that have ALL the reagents
        for reagent_id in reagent_id_list:
            reaction_list = reaction_list.filter(reagents__id=reagent_id)

        reaction_list = reaction_list.select_related('name', 'id', 'process_function', 'solvent').prefetch_related('reagents', 'solvent_properties')
        reaction_data = []
        for reaction in reaction_list:
            reaction_data.append(get_reaction_data(reaction))
    return HttpResponse(json.dumps(reaction_data))

#List links to all reagents
def reagents(request):
    reagent_list = Reagent.objects.all()
    reagent_data = []
    for reagent in reagent_list:
      reagent_data.append(get_reagent_data(reagent))
    return HttpResponse(json.dumps(reagent_data))

#List basic info about a single reagent, id#123 in the database [as JSON]
def reagent(request, id):
    reagent = Reagent.objects.get(id=id)
    reagent_data = get_reagent_data(reagent)
    return HttpResponse(json.dumps(reagent_data))

#If the user entered in [text], what reagent(s) do I get?
def find_reagents(request):
    return HttpResponse(request.GET.get('text', None))

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
#    pseudocode/almost code but it wouldn't actually work
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
