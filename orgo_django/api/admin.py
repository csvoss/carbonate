from django.contrib import admin
from api import models

# Register your models here.
register_models = [models.Property, models.Reagent, models.ReagentSet, models.Reaction]

for model in register_models:
	admin.site.register(model)