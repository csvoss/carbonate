from django.conf.urls import patterns, url

from api import views

SMILESREGEX = '(?P<molecule>.*)'
REACTIONREGEX = '(?P<reaction>.*)'
REAGENTREGEX = '(?P<reagent>.*)'
IDREGEX = '(?P<id>d+)'

urlpatterns = patterns('',
    # ex: /api/
    url(r'^$', views.index, name='index'),
    url(r'^reagents.json/?$', views.reagents_json, name='reagents_json'),
    url(r'^reagents.html/?$', views.reagents_html, name='reagents_html'),
    url(r'^reagents/'+REAGENTREGEX+'/$', views.reagent, name='reagent'),
    url(r'^reactions.json/?$', views.reactions_json, name='reactions_json'),
    url(r'^reactions.html/?$', views.reactions_html, name='reactions_html'),
    url(r'^reactions/'+REACTIONREGEX+'/react/'+SMILESREGEX+'/render/$', views.reaction_reactrender, name='reaction_react_render'),
    url(r'^reactions/'+REACTIONREGEX+'/react/'+SMILESREGEX+'/$', views.reaction_react, name='reaction_react'),
    url(r'^reactions/'+REACTIONREGEX+'/$', views.reaction, name='reaction'),
    url(r'^molecule/'+SMILESREGEX+'/render/$', views.molecule_render, name='molecule_render'),
    url(r'^molecule/'+SMILESREGEX+'/$', views.molecule, name='molecule'),
    url(r'^randomMolecule/render/$', views.random_molecule_render, name='random_molecule_render'),
    url(r'^randomMolecule/$', views.random_molecule, name='random_molecule'),
    url(r'^randomOneStep/$', views.random_one_step, name='random_one_step'),      ## ?rxns = ...
    url(r'^randomSynthesis/$', views.random_synthesis, name='random_synthesis'),    ## ?rxns = ...
    url(r'^synthesisProblems/$', views.synthesis_problems, name='synthesis_problems'),
    url(r'^synthesisProblems/'+IDREGEX+'/reactant/render/$', views.synthesis_problem_reactant_render, name='synthesis_problem_reactant_render'),
    url(r'^synthesisProblems/'+IDREGEX+'/reactant/$', views.synthesis_problem_reactant, name='synthesis_problem_reactant'),
    url(r'^synthesisProblems/'+IDREGEX+'/target/render/$', views.synthesis_problem_target_render, name='synthesis_problem_target_render'),
    url(r'^synthesisProblems/'+IDREGEX+'/target/$', views.synthesis_problem_target, name='synthesis_problem_target'),
    url(r'^synthesisProblems/'+IDREGEX+'/$', views.synthesis_problem, name='synthesis_problem'),
    url(r'^reagentsFromString/$', views.reagents_from_string, name='reagents_from_string'), ## ?str = ...
    url(r'^reactionsxFromReagents/$', views.reactions_from_reagents, name='reactions_from_reagents'),   ## ?reagents = ...
    url(r'^reactionsFromString/$', views.reactions_from_string, name='reactions_from_string'), ## ?str = ...
)

## SMILES testcases
## N1CCN(CC1)C(C(F)=C2)=CC(=C2C4=O)N(C3CC3)C=C4(=O)O