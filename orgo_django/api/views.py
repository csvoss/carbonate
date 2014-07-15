from django.shortcuts import get_object_or_404, render
from django.http import HttpResponseRedirect, HttpResponse
from django.core.urlresolvers import reverse
from engine.renderSVG import render as svg_render
from engine.randomGenerator import randomStart
from engine.toMolecule import moleculify
from engine.toSmiles import smilesify
from engine.helperFunctions import moleculeCompare
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


#List links to all reactions [as JSON]
def reactions(request):
    '''
    Creates a JSON object with all the reactions in url format
    Returns a HTTP Response with the JSON object.
    '''
    reaction_list = Reaction.objects.all().select_related('name', 'id')
    reactions = []
    for reaction in reaction_list:
        name = reaction.name
        r_id = reaction.id
        reactions.append({
            "id": r_id,
            "name": name,
            })
    return HttpResponse(json.dumps(reactions))

#List basic info about a single reaction, id#123 in the database [as JSON]
def reaction(request, id):
    reaction = Reaction.objects.get(id=id)
    reaction_attrs = {
        "id": reaction.id,
        "name": reaction.name,
        "process_function": reaction.process_function,
        "reaction_site_function": reaction.reaction_site_function,
        "reagents": reaction.reagents,
        "solvent": reaction.solvent,
        "solvent_properties": solvent_properties
        }
    return HttpResponse(json.dumps(reaction_attrs))

#If the user entered in [list of reagents, by id], what reaction(s) do I get?
def find_reactions(request):
    if request.method == "GET":
        reagent_id_list = request.GET.get('reagents', None)
        reaction_list = Reaction.objects.all().prefetch_related('reagents__id')
        for reagent_id in reagent_id_list:
            reaction_list = reaction_list.filter(reagents__id=reagent_id)
        reaction_list = reaction_list.select_related('name','id')
        reactions = []
        for reaction in reaction_list:
            name = reaction.name
            r_id = reaction.id
            reactions.append({
                "id": r_id,
                "name": name,
                })
    return HttpResponse(json.dumps(reactions))

#List links to all reagents
def reagents(request):
    reagent_list = Reagent.objects.all()
    reagents = []
    for reagent in reagent_list:
      reagents.append({
        "id":reagent.id,
        "name":reagent.name,
        "is_solvent":reagent.is_solvent,
        "diagram_name":reagent.diagram_name,
        })
    return HttpResponse(json.dumps(reagents))

#List basic info about a single reagent, id#123 in the database [as JSON]
def reagent(request, id):
    reagent = Rereagent.objects.get(id=id)
    attrs = []
    attrs.append({
        "id": id,
        "name": reagent.name,
        "is_solvent": reagent.is_solvent,
        "diagram_name": reagent.diagram_name,
        "smiles": reagent.smiles,
        "properties": reagent.properties,
    })
    return HttpResponse(json.dumps(attrs))

#If the user entered in [text], what reagent(s) do I get?
def find_reagents(request):
    return HttpResponse(request.GET.get('text', None))

#what if the SMILES are different but represent the same molecule?
def check_if_equal (request):
    mol1 = moleculify(request.GET.get('mol1', None))
    mol2 = moleculify(request.GET.get('mol2', None))
    isEqual = moleculeCompare(mol1, mol2)
    return HttpResponse(isEqual)


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
