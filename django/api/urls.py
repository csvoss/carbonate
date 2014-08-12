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
    url(r'^$', views.index, name='index'),
    url(r'^(?i)test/'+SMILESREGEX+r'/$', views.test_parsers, name='test_parsers'),
    url(r'^(?i)reactions/$', views.all_reactions, name='all_reactions'),
    url(r'^(?i)reaction/'+IDREGEX+r'/$', views.get_reaction, name='get_reaction'),
    url(r'^(?i)findReactions/$', views.find_reactions, name='find_reactions'),
    url(r'^(?i)reagents/$', views.all_reagents, name='all_reagents'),
    ## /api/allReagentNames/
    url(r'^(?i)allReagentNames/$', views.all_reagent_names, name='all_reagent_names'),
    url(r'^(?i)reagent/'+IDREGEX+r'/$', views.get_reagent, name='get_reagent'),
    url(r'^(?i)findReagents/$', views.find_reagents, name='find_reagents'),
    url(r'^(?i)checkIfEqual/$', views.check_if_equal, name='check_if_equal'),
    ## /api/isCorrectReagentSet?submitted_reagents=['HBr','HI']&solvent_id=[2]&solution_ids=[2,3]
    url(r'^(?i)isCorrectReagentSet$', views.is_correct_reagent_set, name='is_correct_reagent_set'),
    ## /api/react?reaction=123&reactants=["SM1","SM2"]
    url(r'^(?i)react/$', views.react, name='react'),
    url(r'^(?i)reactSVG/$', views.react_then_SVG, name='react_then_SVG'),
    url(r'^(?i)renderSVG/$', views.render_SVG, name='render_SVG'),
    url(r'^(?i)toCanonical/$', views.to_canonical_view, name='to_canonical_view'),
    url(r'^(?i)randomSmiles/$', views.random_smiles, name='random_smiles'),
    url(r'^(?i)randomSVG/$', views.random_SVG, name='random_SVG'),
)
