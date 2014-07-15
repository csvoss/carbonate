from django.conf.urls import patterns, url

from api import views

SMILESREGEX = r'(?P<smiles>.*)'
IDREGEX = r'(?P<id>.*)'

urlpatterns = patterns(
    '',
    # ex: /api/
    url(r'^$', views.index, name='index'),
    ## /api/test/CCCCCC/test_smiles_to_molecule_and_back
    url(r'^test/'+SMILESREGEX+r'/$', views.test_smiles_to_molecule_and_back, name='test_smiles_to_molecule_and_back'),
    ## /api/reactions
    url(r'^reactions/$', views.reactions, name='reactions'),
    ## /api/reaction/123 
    url(r'^reaction/'+IDREGEX+r'/$', views.reaction, name='reaction'),
    ## /api/findReactions?reagents=[44,2]
    url(r'^findReactions/$', views.find_reactions, name='find_reactions'),
    ## /api/reagents 
    url(r'^reagents/$', views.reagents, name='reagents'),
    ## /api/reagent/123
    url(r'^reagent/'+IDREGEX+r'/$', views.reagent, name='reagent'),
    ## /api/findReagents?text="HBr"
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

## SMILES testcases
## N1CCN(CC1)C(C(F)=C2)=CC(=C2C4=O)N(C3CC3)C=C4(=O)O