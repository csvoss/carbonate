from django.shortcuts import render
from django.shortcuts import get_object_or_404, render
from django.http import HttpResponseRedirect, HttpResponse
from django.core.urlresolvers import reverse
from engine.serverRender import render as svg_render

# Create your views here.

## url(r'^$', views.index, name='index'),
def index(request):
    context = {}
    return render(request, 'api/index.html', context)


## url(r'^reagents.json/?$', views.reagents_json, name='reagents_json'),
def reagents_json(request):
    pass


## url(r'^reagents.html/?$', views.reagents_html, name='reagents_html'),
def reagents_html(request):
    pass


## url(r'^reagents/'+REAGENTREGEX+'/$', views.reagent, name='reagent'),
def reagent(request, reagent):
    pass


## url(r'^reactions.json/?$', views.reactions_json, name='reactions_json'),
def reactions_json(request):
    pass


## url(r'^reactions.html/?$', views.reactions_html, name='reactions_html'),
def reactions_html(request):
    pass


## url(r'^reactions/'+REACTIONREGEX+'/$', views.reaction, name='reaction'),
def reaction(request, reaction):
    pass


## url(r'^reactions/'+REACTIONREGEX+'/react/'+SMILESREGEX+'/$', views.reaction_react, name='reaction_react'),
def reaction_react(request, reaction, molecule):
    pass


## url(r'^reactions/'+REACTIONREGEX+'/react/'+SMILESREGEX+'/render/$', views.reaction_reactrender, name='reaction_react_render'),
def reaction_reactrender(request, reaction, molecule):
    pass


## url(r'^molecule/'+SMILESREGEX+'/$', views.molecule, name='molecule'),
def molecule(request, molecule):
    molecule = molecule.replace('/#','#')
    return HttpResponse(molecule)


## url(r'^molecule/'+SMILESREGEX+'/render/$', views.molecule_render, name='molecule_render'),
def molecule_render(request, molecule):
    molecule = molecule.replace('/#','#')
    return HttpResponse(svg_render(molecule))


## url(r'^randomMolecule/$', views.random_molecule, name='random_molecule'),
def random_molecule(request):
    pass


## url(r'^randomMolecule/render/$', views.random_molecule_render, name='random_molecule_render'),
def random_molecule_render(request):
    pass


## url(r'^randomOneStep/$', views.random_one_step, name='random_one_step'),      ## ?rxns = ...
def random_one_step(request):
    queryparam = request.GET.get('rxns', None)
    pass


## url(r'^randomSynthesis/$', views.random_synthesis, name='random_synthesis'),    ## ?rxns = ...
def random_synthesis(request):
    queryparam = request.GET.get('rxns', None)
    pass


## url(r'^synthesisProblems/$', views.synthesis_problems, name='synthesis_problems'),
def synthesis_problems(request):
    pass


## url(r'^synthesisProblems/'+IDREGEX+'/$', views.synthesis_problem, name='synthesis_problem'),
def synthesis_problem(request):
    pass


## url(r'^synthesisProblems/'+IDREGEX+'/reactant/$', views.synthesis_problem_reactant, name='synthesis_problem_reactant'),
def synthesis_problem_reactant(request):
    pass


## url(r'^synthesisProblems/'+IDREGEX+'/reactant/render/$', views.synthesis_problem_reactant_render, name='synthesis_problem_reactant_render'),
def synthesis_problem_reactant_render(request):
    pass


## url(r'^synthesisProblems/'+IDREGEX+'/target/$', views.synthesis_problem_target, name='synthesis_problem_target'),
def synthesis_problem_target(request):
    pass


## url(r'^synthesisProblems/'+IDREGEX+'/target/render/$', views.synthesis_problem_target_render, name='synthesis_problem_target_render'),
def synthesis_problem_target_render(request):
    pass


## url(r'^reagentsFromString/$', views.reagents_from_string, name='reagents_from_string'), ## ?str = ...
def reagents_from_string(request):
    queryparam = request.GET.get('str', None)
    pass


## url(r'^reactionsxFromReagents/$', views.reactions_from_reagents, name='reactions_from_reagents'),   ## ?reagents = ...
def reactions_from_reagents(request):
    queryparam = request.GET.get('reagents', None)
    pass


## url(r'^reactionsFromString/$', views.reactions_from_string, name='reactions_from_string'), ## ?str = ...
def reactions_from_string(request):
    queryparam = request.GET.get('str', None)
    pass

