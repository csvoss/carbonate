from django.db import models
from api.models import Property, Reagent, ReagentSet, Reaction

class Synthesis(models.Model):
    id = models.AutoField(primary_key=True)
    reactant_smiles = models.TextField(
        help_text = "Starting molecule's SMILES string",
    )
    product_smiles = models.TextField(
        help_text = "Target molecule's SMILES string",
    )
    correct_answer = models.CharField(max_length=100) 

    def __unicode__(self):
        return "Synthesis #%s" % str(self.id)

    class Meta:
        verbose_name_plural = "Syntheses"

class SingleStepProblem(models.Model):
    id = models.AutoField(primary_key=True)
    reactant_smiles = models.TextField(
        help_text = "Reactant's SMILES string",
    )
    product_smiles = models.TextField(
        help_text = "Product's SMILES string",
    )
    correct_answer = models.CharField(max_length=100)

    def __unicode__(self):
        return "SingleStepProblem #%s" % str(self.id)

class SingleStepHardProblem(models.Model):
    id = models.AutoField(primary_key=True)
    reactant_smiles = models.TextField(
        help_text = "Reactant's SMILES string",
    )
    product_smiles = models.TextField(
        help_text = "Product's SMILES string",
    )
    #TODO: change to ManyToManyField to allow multiple solutions
    answer = models.ForeignKey(ReagentSet, help_text="Contains all reagents and solvent for one valid solution")

    def __unicode__(self):
        return "SingleStepHardProblem #%s" % str(self.id)

class PredictProductsProblem(models.Model):
    id = models.AutoField(primary_key=True)
    reactant_smiles = models.TextField(
        help_text = "Reactant's SMILES string",
    )
    reagents = models.CharField(max_length=100, help_text="Human-readable description of reagents added")
    correct_answer = models.TextField(
        help_text="SMILES for the correct product",
    )
    incorrect_answer1 = models.TextField(
        help_text="SMILES for a fake answer",
    )
    incorrect_answer2 = models.TextField(
        help_text="SMILES for a fake answer",
    )
    incorrect_answer3 = models.TextField(
        help_text="SMILES for a fake answer",
    )

    def __unicode__(self):
        return "PredictProductsProblem #%s" % str(self.id)
