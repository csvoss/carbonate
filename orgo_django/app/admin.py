"""
admin.py
Controls which models from app/models.py are accessible from 
the admin interface.
You can access the admin interface at localhost:4000/admin.
"""

from django.contrib import admin
from django.db.models import Model
from app import models
import inspect

## Register all the models!
for i in dir(models):
    thing = getattr(models, i)
    if inspect.isclass(thing):
        if issubclass(thing, Model):
            try:
                admin.site.register(thing)
            except admin.sites.AlreadyRegistered as e:
                pass
