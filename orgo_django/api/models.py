from django.db import models
from fields import StringListField

class Property(models.Model):
    id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=50)
    
class Reagent(models.Model):
    id = models.AutoField(primary_key=True)
    names = StringListField(help_text='All valid names for this reagent')
    isSolvent = models.BooleanField(default=False)
    diagram_name = models.CharField(max_length=50, blank=True, null=True, help_text="HTML-compatible, human-readable name of this reagent")
    smiles = models.CharField(max_length=100, blank=True, null=True,\
        help_text="A SMILES string representation of the molecule")
    properties = models.ManyToManyField(Property, blank=True, null=True, \
        help_text="Useful properties of reagent e.g. aprotic")


class Reaction(models.Model):
    id = models.AutoField(primary_key=True)
    process_function = models.CharField(max_length=50, \
        help_text="Called at each valid place in reactant")
    reaction_site_function = models.CharField(max_length=50, \
        help_text="Finds valid locations in reactant. eg findDoubleBonds")
    reagents = models.ManyToManyField(Reagent, related_name='reactions', help_text="Reagents that trigger this reaction taking place")
    solvent = models.ForeignKey(Reagent, related_name='solventReactions', \
        blank=True, null=True, help_text="Solvents that are compatible with this reaction")
    solvent_properties = models.ManyToManyField(Property, \
        blank=True, null=True)
