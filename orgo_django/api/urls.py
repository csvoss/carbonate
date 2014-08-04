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
    ## /api/reagent/123
    url(r'^reagent/'+IDREGEX+r'/$', views.get_reagent, name='get_reagent'),
    ## /api/getValidReagent?name=hbr
    url(r'^getValidReagent/$', views.get_valid_reagent, name='get_valid_reagent'),
    ## /api/findReagents?name=HBr&properties=["aprotic"]
    url(r'^findReagents/$', views.find_reagents, name='find_reagents'),
    ##/api/checkIfEqual?mol1=ABC&mol2=DEF
    url(r'^checkIfEqual/$', views.check_if_equal, name='check_if_equal'),
    ## /api/react?reaction=123&reactants=["SM1","SM2"]
    url(r'^react/$', views.react, name='react'),
    ## /api/renderSVG?molecule="SMILES"
    url(r'^renderSVG/$', views.render_SVG, name='render_SVG'),
    ## /api/randomGenSmiles/ 
    url(r'^randomGenSmiles/$', views.random_gen_smiles, name='random_gen_smiles'),
    #/api/randomGenSVG
    url(r'^randomGenSVG/$', views.random_gen_SVG, name='random_gen_SVG'),
)
