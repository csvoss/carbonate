from django.conf.urls import patterns, url

from api import views

SMILESREGEX = '(?P<smiles>.*)'
REACTIONREGEX = '(?P<reaction>.*)'
REAGENTREGEX = '(?P<reagent>.*)'
IDREGEX = '(?P<id>d+)'

urlpatterns = patterns('',
    # ex: /api/
    url(r'^$', views.index, name='index'),

    url(r'^molecule/'+SMILESREGEX+r'/render/$', views.molecule_render, name='molecule_render'),
    url(r'^molecule/'+SMILESREGEX+r'/$', views.molecule, name='molecule'),
    
    url(r'^reagents.json/?$', views.reagents_json, name='reagents_json'),
    url(r'^reagents.html/?$', views.reagents_html, name='reagents_html'),
    url(r'^reagents/'+REAGENTREGEX+r'/$', views.reagent, name='reagent'),
    
    url(r'^reactions.json/?$', views.reactions_json, name='reactions_json'),
    url(r'^reactions.html/?$', views.reactions_html, name='reactions_html'),
    url(r'^reactions/'+REACTIONREGEX+r'/react/'+SMILESREGEX+r'/render/$', views.reaction_reactrender, name='reaction_react_render'),
    url(r'^reactions/'+REACTIONREGEX+r'/react/'+SMILESREGEX+r'/$', views.reaction_react, name='reaction_react'),
    url(r'^reactions/'+REACTIONREGEX+r'/$', views.reaction, name='reaction'),

    url(r'^randomMolecule/render/$', views.random_molecule_render, name='random_molecule_render'),
    url(r'^randomMolecule/$', views.random_molecule, name='random_molecule'),
    url(r'^randomOneStep/$', views.random_one_step, name='random_one_step'),      ## ?rxns = ...
    url(r'^randomSynthesis/$', views.random_synthesis, name='random_synthesis'),    ## ?rxns = ...

    url(r'^synthesisProblems/$', views.synthesis_problems, name='synthesis_problems'),
    url(r'^synthesisProblems/'+IDREGEX+r'/reactant/render/$', views.synthesis_problem_reactant_render, name='synthesis_problem_reactant_render'),
    url(r'^synthesisProblems/'+IDREGEX+r'/reactant/$', views.synthesis_problem_reactant, name='synthesis_problem_reactant'),
    url(r'^synthesisProblems/'+IDREGEX+r'/target/render/$', views.synthesis_problem_target_render, name='synthesis_problem_target_render'),
    url(r'^synthesisProblems/'+IDREGEX+r'/target/$', views.synthesis_problem_target, name='synthesis_problem_target'),
    url(r'^synthesisProblems/'+IDREGEX+r'/$', views.synthesis_problem, name='synthesis_problem'),

    url(r'^reagentsFromString/$', views.reagents_from_string, name='reagents_from_string'), ## ?str = ...
    url(r'^reactionsxFromReagents/$', views.reactions_from_reagents, name='reactions_from_reagents'),   ## ?reagents = ...
    url(r'^reactionsFromString/$', views.reactions_from_string, name='reactions_from_string'), ## ?str = ...


    url(r'^test/'+SMILESREGEX+r'/$', views.test_smiles_to_molecule_and_back, name='test_smiles_to_molecule_and_back'),
)

## SMILES testcases
## N1CCN(CC1)C(C(F)=C2)=CC(=C2C4=O)N(C3CC3)C=C4(=O)O