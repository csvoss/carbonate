from django.db import models
from fields import StringListField

class Reagent(models.Model):
	id = models.AutoField(primary_key=True)
	names = StringListField(help_text='All valid names for this reagent')
	isSolvent = models.BooleanField(default=False)
	diagram_name = models.CharField(max_length=50, blank=True, null=True)
	smiles = models.CharField(max_length=100, blank=True, null=True,\
		help_text="A string representation of the molecule")
	properties = models.ManyToManyField(Property, blank=True, null=True, \
		help_text="Useful properties of reagent e.g. aprotic")


class Reaction(models.Model):
	id = models.AutoField(primary_key=True)
	process_function = models.CharField(max_length=50, \
		help_text="Called at each valid place in reactant")
	reaction_site_function = models.CharField(max_length=50, \
		help_text="Finds valid locations in reactant. eg findDoubleBonds")
	reagents = models.ManyToManyField(Reagent, related_name='reactions')
	solvent = models.ForeignKey(Reagent, related_name='solventReactions', \
		blank=True, null=True)
	solvent_properties = models.ManyToManyField(Property, \
		blank=True, null=True)

class Property(models.Model):
	id = models.AutoField(primary_key=True)
	name = models.CharField(max_length=50)