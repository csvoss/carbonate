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
    url(r'^test/'+SMILESREGEX+r'/$', views.test_parsers, name='test_parsers'),
    url(r'^reactions/$', views.all_reactions, name='all_reactions'),
    url(r'^reaction/'+IDREGEX+r'/$', views.get_reaction, name='get_reaction'),
    url(r'^findReactions/$', views.find_reactions, name='find_reactions'),
    url(r'^reagents/$', views.all_reagents, name='all_reagents'),
    ## /api/allReagentNames/
    url(r'^allReagentNames/$', views.all_reagent_names, name='all_reagent_names'),
    url(r'^reagent/'+IDREGEX+r'/$', views.get_reagent, name='get_reagent'),
    url(r'^findReagents/$', views.find_reagents, name='find_reagents'),
    url(r'^checkIfEqual/$', views.check_if_equal, name='check_if_equal'),
    ## /api/isCorrectReagentSet?submitted_reagents=['HBr','HI']&solvent_id=[2]&solution_ids=[2,3]
    url(r'^isCorrectReagentSet$', views.is_correct_reagent_set, name='is_correct_reagent_set'),
    ## /api/react?reaction=123&reactants=["SM1","SM2"]
    url(r'^react/$', views.react, name='react'),
    url(r'^reactSVG/$', views.react_then_SVG, name='react_then_SVG'),
    url(r'^renderSVG/$', views.render_SVG, name='render_SVG'),
    url(r'^toCanonical/$', views.to_canonical_view, name='to_canonical_view'),
    url(r'^randomSmiles/$', views.random_smiles, name='random_smiles'),
    url(r'^randomSVG/$', views.random_SVG, name='random_SVG'),
    url(r'^updateReagentNames/$', views.update_reagent_names, name='update_reagent_names'),
)
