from .base import *

from mongoengine import connect

# Database
# https://docs.djangoproject.com/en/1.9/ref/settings/#databases

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': 'aagdb',
        'USER': 'aagdb',
        'PASSWORD': 'No5waip6guw_e2ie',
        'HOST': 'db-user.rc.fas.harvard.edu',
        'PORT': '5432',
    },

    'mongo': {
       'ENGINE': 'django.db.backends.dummy',
       'NAME': 'test',
       'USER': 'test',
       'PASSWORD': 'test',
       'HOST': 'molspace.rc.fas.harvard.edu',
       'PORT': '',
       'ALIAS': 'default',
       'PROJECTS': ['cep']
    }
}

mongo = DATABASES["mongo"]
connect(mongo["NAME"],
        alias=mongo['ALIAS'],
        host=mongo["HOST"],
        username=mongo["USER"],
        password=mongo["PASSWORD"])
