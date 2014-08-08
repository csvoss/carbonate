"""
api/urls

URLs for the API views.
"""

from django.conf.urls import patterns, url

from api import views

SMILESREGEX = r'(?P<smiles>.*)'
IDREGEX = r'(?P<pk>.*)'

urlpatterns = patterns(
    '',
    # ex: /api/
    url(r'^$', views.index, name='index'),
    ## /api/test/CCCCCC/
    url(r'^test/'+SMILESREGEX+r'/$', views.test_parsers, name='test_parsers'),
    ## /api/reactions
    url(r'^reactions/$', views.all_reactions, name='all_reactions'),
    ## /api/reaction/123 
    url(r'^reaction/'+IDREGEX+r'/$', views.get_reaction, name='get_reaction'),
    ## /api/findReactions?reagents=[44,2]
    url(r'^findReactions/$', views.find_reactions, name='find_reactions'),
    ## /api/reagents 
    url(r'^reagents/$', views.all_reagents, name='all_reagents'),
    ## /api/allReagentNames/
    url(r'^allReagentNames/$', views.all_reagent_names, name='all_reagent_names'),
    ## /api/reagent/123
    url(r'^reagent/'+IDREGEX+r'/$', views.get_reagent, name='get_reagent'),
    ## /api/findReagents?name=HBr&properties=["aprotic"]
    url(r'^findReagents/$', views.find_reagents, name='find_reagents'),
    ##/api/checkIfEqual?mol1=ABC&mol2=DEF
    url(r'^checkIfEqual/$', views.check_if_equal, name='check_if_equal'),
    ## /api/isCorrectReagentSet?submitted_reagents=['HBr','HI']&solvent_id=[2]&solution_ids=[2,3]
    url(r'^isCorrectReagentSet$', views.is_correct_reagent_set, name='is_correct_reagent_set'),
    ## /api/react?reaction=123&reactants=["SM1","SM2"]
    url(r'^react/$', views.react, name='react'),
    ## /api/renderSVG?molecule="SMILES"
    url(r'^renderSVG/$', views.render_SVG, name='render_SVG'),
    ## /api/toCanonical?molecule="SMILES"
    url(r'^toCanonical/$', views.to_canonical_view, name='to_canonical_view'),
    ## /api/randomGenSmiles/ 
    url(r'^randomGenSmiles/$', views.random_gen_smiles, name='random_gen_smiles'),
    #/api/randomGenSVG
    url(r'^randomGenSVG/$', views.random_gen_SVG, name='random_gen_SVG'),
)
