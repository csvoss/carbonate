from django.conf.urls import patterns, url
from app import views

urlpatterns = patterns('',
    # ex: /api/
    url(r'^$', views.index, name='index'),
)

