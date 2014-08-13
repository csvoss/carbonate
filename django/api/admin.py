"""
admin.py
Controls which models from api/models.py are accessible from 
the admin interface.
You can access the admin interface at localhost:4000/admin.
"""

from django.contrib import admin
from api import models

def register(model):
    def make_class(admin_class):
        admin.site.register(model, admin_class)
        return admin_class
    return make_class


class ReactionForReagentSetInline(admin.TabularInline):
    model = models.Reaction
    fk_name = 'reagent_set'
    verbose_name = "Reaction"
    verbose_name_plural = "Reactions which use this reagent set"

class ReagentSetForSolventInline(admin.TabularInline):
    model = models.ReagentSet
    fk_name = 'solvent'
    verbose_name = "Reagent set"
    verbose_name_plural = "Reagent sets which use this as a solvent"
    raw_id_fields = ('reagents', 'solvent_properties')

class PropertyForReagentInline(admin.TabularInline):
    model = models.Reagent.properties.through
    verbose_name = "Property"
    verbose_name_plural = "Properties of this reagent"

class SolventPropertiesForReagentSetInline(admin.TabularInline):
    model = models.ReagentSet.solvent_properties.through
    verbose_name = "Solvent property"
    verbose_name_plural = "Solvent properties of this reagent set"

class ReagentsForReagentSetInline(admin.TabularInline):
    model = models.ReagentSet.reagents.through
    verbose_name = "Reagent"
    verbose_name_plural = "Reagents of this reagent set"

@register(models.Reagent)
class ReagentAdmin(admin.ModelAdmin):
    fieldsets = (
        ('Options', {
            'classes': ('wide', 'extrapretty'),
            'fields': ('diagram_name', 'is_solvent', 'smiles', 'properties')
        }),
    )
    def to_string(self, obj):
        return str(obj)
    list_display = ('to_string', 'diagram_name', 'is_solvent', 'smiles')
    list_display_links = ('to_string',)
    list_editable = ('diagram_name', 'is_solvent', 'smiles')
    list_filter = ('is_solvent',)
    search_fields = ['diagram_name', 'smiles']
    raw_id_fields = ('properties',)
    view_on_site = True
    # inlines = [ReagentSetForSolventInline, PropertyForReagentInline]

@register(models.ReagentName)
class ReagentNameAdmin(admin.ModelAdmin):
    pass

@register(models.Property)
class PropertyAdmin(admin.ModelAdmin):
    fieldsets = (
        ('Options', {
            'classes': ('wide', 'extrapretty'),
            'fields': ('name',)
        }),
    )
    def to_string(self, obj):
        return str(obj)
    list_display = ('to_string', 'name', 'id')
    list_display_links = ('to_string',)
    list_editable = ('name', 'id')
    search_fields = ['name']

@register(models.ReagentSet)
class ReagentSetAdmin(admin.ModelAdmin):
    fieldsets = (
        ('Options', {
            'classes': ('wide', 'extrapretty'),
            'fields': ('name', 'reagents', 'solvent', 'solvent_properties')
        }),
    )
    def to_string(self, obj):
        return str(obj)
    to_string.short_description = 'Reagent set'
    list_display = ('to_string', 'name', 'solvent', 'id')
    list_display_links = ('to_string',)
    list_editable = ('name', 'solvent',)
    list_filter = ('reagents', 'solvent')
    raw_id_fields = ('reagents', 'solvent_properties')
    # formfield_overrides = {
    #     ManyToManyField: {'widget': CheckboxSelectMultiple},  ## checkboxes
    # }
    search_fields = ['name', 'reagents__name', 'solvent__name']

    # inlines = [
    #     ReactionForReagentSetInline,
    #     ReagentsForReagentSetInline,
    #     SolventPropertiesForReagentSetInline,
    # ]


@register(models.Reaction)
class ReactionAdmin(admin.ModelAdmin):
    fieldsets = (
        ('Options', {
            'classes': ('wide', 'extrapretty'),
            'fields': ('name', 'process_function', 'reagent_set')
        }),
    )
    def to_string(self, obj):
        return str(obj)

    list_display = ('to_string', 'name', 'id')
    list_display_links = ('to_string',)
    radio_fields = {'reagent_set': admin.VERTICAL}
    view_on_site = True

