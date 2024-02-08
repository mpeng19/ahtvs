from .base import *

from mongoengine import connect

DEBUG = True

# Database
# https://docs.djangoproject.com/en/1.9/ref/settings/#databases

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': 'postgres',
        'USER': 'postgres',
        'PASSWORD': 'postgres',
        'HOST': 'localhost',
        'PORT': '5432',
    },

    'mongo': {
       'ENGINE': 'django.db.backends.dummy',
       'NAME': 'flow_batt',
       'USER': 'flow_batt',
       'PASSWORD': 'muchpotential!',
       'HOST': 'molspace.rc.fas.harvard.edu',
       'PORT': '',
       'ALIAS': 'default',
       'PROJECTS': ['flow_batt']
    }
}


mongo = DATABASES["mongo"]
connect(mongo["NAME"], alias=mongo["ALIAS"], host=mongo["HOST"], username=mongo["USER"], password=mongo["PASSWORD"])
