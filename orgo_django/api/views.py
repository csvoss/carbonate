from django.shortcuts import render
from django.shortcuts import get_object_or_404, render
from django.http import HttpResponseRedirect, HttpResponse
from django.core.urlresolvers import reverse
from engine.renderSVG import render as svg_render
from engine.toMolecule import moleculify
from engine.toSmiles import smilesify

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

def reactions(request):
    return HttpResponse("hi")

def reaction(request, id):
    return HttpResponse("hi")

def find_reactions(request):
    return HttpResponse(request.GET.get('reagents', None))

def reagents(request):
    return HttpResponse("hi")

def reagent(request, id):
    agent = Reagent.objects.get(id=id)
    return HttpResponse(json(test_smiles_to_molecule_and_back(request,agent)))

def find_reagents(request):
    return HttpResponse(request.GET.get('text', None))

def react(request):
    return HttpResponse(request.GET.get('reaction', None) + '\n' + request.GET.get('molecules', None))

def render_SVG(request):
    smiles = request.GET.get('molecule', None) # change to error raise later
    hydrogens = request.GET.get('hydrogens', False)
    return HttpResponse(svg_render(smiles, hydrogens)) # hydrogens default to false

def random_gen_smiles(request):
    return HttpResponse("hi")
