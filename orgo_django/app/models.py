from django.db import models
from api.fields import StringListField

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

class PredictProductsProblem(models.Model):
    id = models.AutoField(primary_key=True)
    reactantSmiles = models.CharField(max_length=100)
    reagents = models.CharField(max_length=100)
    correctAnswer = models.CharField(max_length=100)
    incorrectAnswer1 = models.CharField(max_length=100)
    incorrectAnswer2 = models.CharField(max_length=100)
    incorrectAnswer3 = models.CharField(max_length=100)
