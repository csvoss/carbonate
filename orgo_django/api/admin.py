"""
admin.py
Controls which models from api/models.py are accessible from 
the admin interface.
You can access the admin interface at localhost:4000/admin.
"""

from django.contrib import admin
from django.db.models import Model
from api import models
import inspect

## Register ALL the models!
# for i in dir(models):
#     thing = getattr(models, i)
#     if inspect.isclass(thing):
#         if issubclass(thing, Model):
#             try:
#                 admin.site.register(thing)
#             except admin.sites.AlreadyRegistered as e:
#                 pass

def register(model):
    def make_class(admin_class):
        admin.site.register(model, admin_class)
        return admin_class
    return make_class

@register(models.Reagent)
class ReagentAdmin(admin.ModelAdmin):
    fieldsets = (
        (None, {
            'classes': ('wide', 'extrapretty'),
            'fields': ('name','diagram_name')
        }),
        ('Advanced options', {
            'classes': ('wide', 'extrapretty', 
                'collapse'),
            'fields': ('is_solvent', 'smiles', 'properties')
        }),
    )
    list_display = ('name', 'diagram_name', 'is_solvent', 'smiles')

@register(models.Property)
class PropertyAdmin(admin.ModelAdmin):
    pass

@register(models.ReagentSet)
class ReagentSetAdmin(admin.ModelAdmin):
    pass

@register(models.Reaction)
class ReactionAdmin(admin.ModelAdmin):
    pass


