from django.shortcuts import get_object_or_404, render
from django.http import HttpResponseRedirect, HttpResponse
from django.core.urlresolvers import reverse
from engine.renderSVG import render as svg_render
from engine.randomGenerator import randomStart
from engine.toMolecule import moleculify
from engine.toSmiles import smilesify

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
#
#/api/findReagents?text="HBr"
#
#/api/react?reaction=123&molecules=["SM1","SM2"]
#/api/renderSVG?molecule="SMILES"

#List links to all reactions [as JSON]
def reactions(request):
    '''
    Creates a JSON object with all the reactions in url format
    Returns a HTTP Response with the JSON object.
    '''
    reaction_list = Reaction.objects.all()
    reactions = []
    for reaction in reaction_list:
        name = reaction.name
        r_id = reaction.id
        reactions.append({
            "id": r_id,
            "name": name,
            "url":'reaction/'+r_id+'/',
            })
    return HttpResponse(json.dumps(reactions))

#List basic info about a single reaction, id#123 in the database [as JSON]
def reaction(request, id):
    return HttpResponse("hi")

#If the user entered in [list of reagents, by id], what reaction(s) do I get?
def find_reactions(request):
    return HttpResponse(request.GET.get('reagents', None))

#List links to all reagents
def reagents(request):
    return HttpResponse("hi")

#List basic info about a single reagent, id#123 in the database [as JSON]
def reagent(request, id):
    agent = Reagent.objects.get(id=id)
    attrs = []
    attrs.append({"id":id,"name":agent.name,"is_solvent":agent.is_solvent,"diagram_name",agent.diagram_name,"smiles":agent.smiles,"properties":agent.properties})
    return HttpResponse(json.dumps(attrs))

#If the user entered in [text], what reagent(s) do I get?
def find_reagents(request):
    return HttpResponse(request.GET.get('text', None))

#React molecule(s) (by SMILES) with a particular reaction, and return the result (as SMILES)
def react(request):
#    pseudocode/almost code but it wouldn't actually work
#    reactionID = request.GET.get('reaction', None)
#    reactants = request.GET.get('reactants', None)
#    reaction = getFromDatabaseReactionByID(reactionID)
#    intermediate = reactants
#    for i in reaction.locationOfExecution(reactants):
#        intermediate = reaction.function(intermediate)
#    product = intermediate
#    return HttpResponse(product)
    return HttpResponse(request.GET.get('reaction', None) + '\n' + request.GET.get('reactants', None))

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
