from django.conf.urls import patterns, url

from api import views

SMILESREGEX = r'(?P<smiles>.*)'
REACTIONREGEX = r'(?P<reaction>.*)'
REAGENTREGEX = r'(?P<reagent>.*)'
IDREGEX = r'(?P<id>d+)'

urlpatterns = patterns('',
    # ex: /api/
    url(r'^$', views.index, name='index'),

    url(r'^molecule/'+SMILESREGEX+r'/render/$', views.molecule_render, name='molecule_render'),
    url(r'^test/'+SMILESREGEX+r'/$', views.test_smiles_to_molecule_and_back, name='test_smiles_to_molecule_and_back'),
)

## SMILES testcases
## N1CCN(CC1)C(C(F)=C2)=CC(=C2C4=O)N(C3CC3)C=C4(=O)O