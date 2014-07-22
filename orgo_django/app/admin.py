from django.contrib import admin
from app import models

# Register your models here.
register_models = [models.SingleStepProblem, models.SingleStepHardProblem, models.PredictProductsProblem]

for model in register_models:
	admin.site.register(model)
