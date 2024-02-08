from django.contrib import admin
from django.apps import apps

pgmols = apps.get_app_config('pgmols')
for model in pgmols.get_models():
    admin.site.register(model)
