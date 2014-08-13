from django.conf.urls import patterns, url
from app import views

IDREGEX = r'(?P<id>\d+)'

urlpatterns = patterns(
    '',
    url(r'^$', views.index, name='index'),
    url(r'^(?i)synthesis/'+IDREGEX+r'/$', views.synthesis, name='synthesis'),
    url(r'^(?i)singleStep/'+IDREGEX+r'/$', views.single_step, name='single_step'),
    url(r'^(?i)singleStepHard/'+IDREGEX+r'/$', views.single_step_hard, name='single_step'),
    url(r'^(?i)predictProducts/'+IDREGEX+r'/$', views.predict_products, name='predict_products'),
    url(r'^(?i)reactionTutor/'+IDREGEX+r'/$', views.reaction_tutorial, name='reaction_tutorial'),
)
