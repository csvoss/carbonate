from django.db import models
from api.fields import StringListField

class SingleStepProblem(models.Model):
    
    id = models.AutoField(primary_key=True)
    reactantSmiles = models.CharField(max_length=100)
    productSmiles = models.CharField(max_length=100)
    correctAnswer = models.CharField(max_length=100)
    