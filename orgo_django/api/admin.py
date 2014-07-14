from django.contrib import admin
from api import models

# Register your models here.
admin.site.register(models.Property)
admin.site.register(models.Reagent)
admin.site.register(models.Reaction)