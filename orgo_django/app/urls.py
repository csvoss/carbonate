from django.conf.urls import patterns, url
from app import views

urlpatterns = patterns(
    '',
    # ex: /api/
    url(r'^$', views.index, name='index'),
    
#    app/synthesis
    url(r'^synthesis/$', views.synthesis, name='synthesis'),
#    app/singleStep
    url(r'^singleStep/$', views.single_step, name='single_step'),
#    app/predictProducts
    url(r'^predictProducts/$', views.predict_products, name='predict_products'),
                       
)

