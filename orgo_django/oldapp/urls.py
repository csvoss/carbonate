from django.conf.urls import patterns, include, url 
from django.contrib.staticfiles.urls import staticfiles_urlpatterns
from django.conf import settings

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('',
    # Examples:
    url(r'^[oldapp/]*$', 'oldapp.views.home', name='home'),
    url(r'^[oldapp/]*api/signup/$', 'oldapp.views.signUp', name='signUp'),
#    url(r'^[oldapp/]*api/home/$', 'oldapp.views.returnToLoggedInHome', name='returnToLoggedInHome'),
    url(r'^[oldapp/]*api/login/$', 'oldapp.views.logIn', name='logIn'),
    url(r'^[oldapp/]*api/returnReagentHtml/$', 'oldapp.views.makeReagentHtml', name='makeReagentHtml'),
    url(r'^[oldapp/]*api/resetPW/$', 'oldapp.views.resetPW', name='resetPW'),
    url(r'^[oldapp/]*api/changePW/$', 'oldapp.views.changePW', name='changePW'),
    url(r'^[oldapp/]*logout/$', 'oldapp.views.logOut', name='logOut'),
    url(r'^[oldapp/]*api/problemInterface/$', 'oldapp.views.renderProblem', name='renderProblem'),    
    url(r'^[oldapp/]*api/checkSingleStepReaction/$', 'oldapp.views.checkNameReagent', name='checkNameReagent'), 
    url(r'^[oldapp/]*api/showSingleStepAnswer/$', 'oldapp.views.showNRAnswer', name='showNRAnswer'),   
    url(r'^[oldapp/]*api/outpsmiles/$', 'oldapp.views.outpSmiles', name='outpSmiles'),  ###Can delete; this is me learning Django  
    url(r'^[oldapp/]*homeMoleculeChanger/$', 'oldapp.views.homeMoleculeChanger', name='homeMoleculeChanger'),
    url(r'^[oldapp/]*chat/helperpoll/$', 'oldapp.views.helperWaitPoll', name='helperWaitPoll'),
    url(r'^[oldapp/]*chat/helpeepoll/$', 'oldapp.views.helpeeWaitPoll', name='helpeeWaitPoll'),
    url(r'^[oldapp/]*chat/askforhelp/$', 'oldapp.views.askForHelp', name='askForHelp'),
    url(r'^[oldapp/]*chat/volunteertohelp/$', 'oldapp.views.volunteerToHelp', name='volunteerToHelp'),
    url(r'^[oldapp/]*chat/helpeechatpoll/$', 'oldapp.views.helpeeChatPoll', name='helpeeChatPoll'),
    url(r'^[oldapp/]*chat/helperchatpoll/$', 'oldapp.views.helperChatPoll', name='helperChatPoll'),
    
    url(r'^[oldapp/]*renderProblem/$', 'oldapp.views.renderProblem', name='renderProblem'),
    
    url(r'^[oldapp/]*renderSynthesis/$', 'oldapp.views.renderSynthesis', name='renderSynthesis'),
    url(r'^[oldapp/]*renderSynthesis/resume/$', 'oldapp.views.renderOldSynthesis', name='renderOldSynthesis'),
    url(r'^[oldapp/]*namereagent/$', 'oldapp.views.renderNameReagent', name='renderNameReagent'),
    url(r'^[oldapp/]*namereagent/resume/$', 'oldapp.views.renderOldNameReagent', name='renderOldNameReagent'),
                       
                       
    url(r'^[oldapp/]*reactions/$', 'oldapp.views.renderReactions', name='renderReactions'),
    url(r'^[oldapp/]*faq/$', 'oldapp.views.renderFaq', name='renderFaq'),
    
    
    url(r'^[oldapp/]*api/getSynthesisData/$', 'oldapp.views.getSynthesisData', name = 'getSynthesisData'),
    url(r'^[oldapp/]*api/displaySolution/$', 'oldapp.views.getSolutionData', name = 'getSolutionData'),
    url(r'^[oldapp/]*api/addMoleculeToMolecule/$', 'oldapp.views.addMoleculeToMolecule', name = 'addMoleculeToMolecule'),
    url(r'^[oldapp/]*api/addReagentToMolecule/$', 'oldapp.views.addReagentToMolecule', name = 'addReagentToMolecule'),
    url(r'^[oldapp/]*api/deleteMolecule/$', 'oldapp.views.deleteMolecule', name = 'deleteMolecule'),
    url(r'^[oldapp/]*api/saveProblem/$', 'oldapp.views.saveProblem', name = 'saveProblem'),
    url(r'^[oldapp/]*loadSynthesisFromId/$', 'oldapp.views.loadSynthesisFromId', name = 'loadSynthesisFromId'),
    

    # Uncomment the admin/doc line below to enable admin documentation:
    url(r'^[oldapp/]*admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    url(r'^[oldapp/]*admin/', include(admin.site.urls)),
)





if settings.DEBUG:
    urlpatterns += patterns('',
        (r'^[oldapp/]*static/(?P<path>.*)$', 'django.views.static.serve', {'document_root': settings.STATIC_ROOT}),
    )

