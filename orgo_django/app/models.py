from django.db import models
from api.models import Property, Reagent, ReagentSet, Reaction

class Synthesis(models.Model):
    id = models.AutoField(primary_key=True)
    reactantSmiles = models.CharField(max_length=100)
    productSmiles = models.CharField(max_length=100)
    correctAnswer = models.CharField(max_length=100)

class SingleStepProblem(models.Model):
    id = models.AutoField(primary_key=True)
    reactantSmiles = models.CharField(max_length=100)
    productSmiles = models.CharField(max_length=100)
    correctAnswer = models.CharField(max_length=100)

class SingleStepHardProblem(models.Model):
    id = models.AutoField(primary_key=True)
    reactant_smiles = models.TextField()
    product_smiles = models.TextField()
    #TODO: change to ManyToManyField to allow multiple solutions
    answer = models.ForeignKey(ReagentSet, help_text="Contains all reagents and solvent for one valid solution")

    def __unicode__(self):
        return "SingleStepHardProblem #"+str(self.id)

class PredictProductsProblem(models.Model):
    id = models.AutoField(primary_key=True)
    reactantSmiles = models.CharField(max_length=100)
    reagents = models.CharField(max_length=100)
    correctAnswer = models.CharField(max_length=100)
    incorrectAnswer1 = models.CharField(max_length=100)
    incorrectAnswer2 = models.CharField(max_length=100)
    incorrectAnswer3 = models.CharField(max_length=100)
