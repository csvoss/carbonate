from django.conf.urls import patterns, include, url 
from django.contrib.staticfiles.urls import staticfiles_urlpatterns
from django.conf import settings

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('',
    # Examples:
    url(r'^[orgo/]*$', 'orgo.views.home', name='home'),
    url(r'^[orgo/]*api/signup/$', 'orgo.views.signUp', name='signUp'),
#    url(r'^[orgo/]*api/home/$', 'orgo.views.returnToLoggedInHome', name='returnToLoggedInHome'),
    url(r'^[orgo/]*api/login/$', 'orgo.views.logIn', name='logIn'),
    url(r'^[orgo/]*api/returnReagentHtml/$', 'orgo.views.makeReagentHtml', name='makeReagentHtml'),
    url(r'^[orgo/]*api/resetPW/$', 'orgo.views.resetPW', name='resetPW'),
    url(r'^[orgo/]*api/changePW/$', 'orgo.views.changePW', name='changePW'),
    url(r'^[orgo/]*logout/$', 'orgo.views.logOut', name='logOut'),
    url(r'^[orgo/]*api/problemInterface/$', 'orgo.views.renderProblem', name='renderProblem'),    
    url(r'^[orgo/]*api/checkSingleStepReaction/$', 'orgo.views.checkNameReagent', name='checkNameReagent'), 
    url(r'^[orgo/]*api/showSingleStepAnswer/$', 'orgo.views.showNRAnswer', name='showNRAnswer'),   
    url(r'^[orgo/]*api/outpsmiles/$', 'orgo.views.outpSmiles', name='outpSmiles'),  ###Can delete; this is me learning Django  
    url(r'^[orgo/]*homeMoleculeChanger/$', 'orgo.views.homeMoleculeChanger', name='homeMoleculeChanger'),
    url(r'^[orgo/]*chat/helperpoll/$', 'orgo.views.helperWaitPoll', name='helperWaitPoll'),
    url(r'^[orgo/]*chat/helpeepoll/$', 'orgo.views.helpeeWaitPoll', name='helpeeWaitPoll'),
    url(r'^[orgo/]*chat/askforhelp/$', 'orgo.views.askForHelp', name='askForHelp'),
    url(r'^[orgo/]*chat/volunteertohelp/$', 'orgo.views.volunteerToHelp', name='volunteerToHelp'),
    url(r'^[orgo/]*chat/helpeechatpoll/$', 'orgo.views.helpeeChatPoll', name='helpeeChatPoll'),
    url(r'^[orgo/]*chat/helperchatpoll/$', 'orgo.views.helperChatPoll', name='helperChatPoll'),
    
    url(r'^[orgo/]*renderProblem/$', 'orgo.views.renderProblem', name='renderProblem'),
    
    url(r'^[orgo/]*renderSynthesis/$', 'orgo.views.renderSynthesis', name='renderSynthesis'),
    url(r'^[orgo/]*renderSynthesis/resume/$', 'orgo.views.renderOldSynthesis', name='renderOldSynthesis'),
    url(r'^[orgo/]*namereagent/$', 'orgo.views.renderNameReagent', name='renderNameReagent'),
    url(r'^[orgo/]*namereagent/resume/$', 'orgo.views.renderOldNameReagent', name='renderOldNameReagent'),
                       
                       
    url(r'^[orgo/]*reactions/$', 'orgo.views.renderReactions', name='renderReactions'),
    url(r'^[orgo/]*faq/$', 'orgo.views.renderFaq', name='renderFaq'),
    
    
    url(r'^[orgo/]*api/getSynthesisData/$', 'orgo.views.getSynthesisData', name = 'getSynthesisData'),
    url(r'^[orgo/]*api/displaySolution/$', 'orgo.views.getSolutionData', name = 'getSolutionData'),
    url(r'^[orgo/]*api/addMoleculeToMolecule/$', 'orgo.views.addMoleculeToMolecule', name = 'addMoleculeToMolecule'),
    url(r'^[orgo/]*api/addReagentToMolecule/$', 'orgo.views.addReagentToMolecule', name = 'addReagentToMolecule'),
    url(r'^[orgo/]*api/deleteMolecule/$', 'orgo.views.deleteMolecule', name = 'deleteMolecule'),
    url(r'^[orgo/]*api/saveProblem/$', 'orgo.views.saveProblem', name = 'saveProblem'),
    url(r'^[orgo/]*loadSynthesisFromId/$', 'orgo.views.loadSynthesisFromId', name = 'loadSynthesisFromId'),
    

    # Uncomment the admin/doc line below to enable admin documentation:
    url(r'^[orgo/]*admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    url(r'^[orgo/]*admin/', include(admin.site.urls)),
)





if settings.DEBUG:
    urlpatterns += patterns('',
        (r'^[orgo/]*static/(?P<path>.*)$', 'django.views.static.serve', {'document_root': settings.STATIC_ROOT}),
    )

