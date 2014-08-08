from django.db import models
from api.fields import StringListField
from api.engine.toCanonical import to_canonical

class Property(models.Model):
    id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=50)

    class Meta:
        verbose_name_plural = "Properties"

    def __unicode__(self):
        return self.name


class Reagent(models.Model):
    id = models.AutoField(primary_key=True)
    # names = StringListField(help_text='All valid names for this reagent')
    is_solvent = models.BooleanField(default=False)
    diagram_name = models.CharField(max_length=50, blank=True, help_text="HTML-compatible, human-readable name of this reagent")
    smiles = models.CharField(max_length=100, help_text="A SMILES string representation of the molecule (optional)")
    properties = models.ManyToManyField(Property, blank=True, null=True, help_text="Useful properties of reagent e.g. aprotic")

    def __unicode__(self):
        return to_canonical(self.smiles)

class ReagentName(models.Model):
    id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=100, unique=True)
    reagent = models.ForeignKey(Reagent, related_name="names", help_text="The reagent that this name belongs to.")

    def __unicode__(self):
        return self.name

class ReagentSet(models.Model):
    id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=50, blank=True, help_text="Optional name for this set of reagents.")
    reagents = models.ManyToManyField(Reagent, related_name="reagentSets")
    solvent = models.ForeignKey(Reagent, blank=True, null=True, help_text="Use if only specific solvents are used with this reaction set")
    solvent_properties = models.ManyToManyField(Property, blank=True, null=True, help_text="Use if any reagent with given properties works with this set.")

    def __unicode__(self):
        if len(self.name) > 1:
            return self.name
        reagent_names = ""
        for reagent in self.reagents.all():
            reagent_names += str(reagent.name) + ", "
        reagent_names = reagent_names.rstrip(', ')
        solvent_names = ""
        if self.solvent is not None:
            solvent_names += str(self.solvent.name)
        if self.solvent_properties is not None:
            for prop in self.solvent_properties.all().select_related("name"):
                solvent_names += str(prop.name) + ", "
        solvent_names = solvent_names.rstrip(', ')

        display = reagent_names + " Solvent: " + solvent_names
        return unicode(display)


class Reaction(models.Model):
    id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=100)
    process_function = models.CharField(max_length=50, help_text="Called at each valid place in reactant")
    reagent_set = models.ForeignKey(ReagentSet, related_name='reactions', help_text="Reagent set that triggers this reaction")

    def __unicode__(self):
        return self.name

