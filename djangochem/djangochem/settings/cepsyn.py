from .base import *

from mongoengine import connect

# Database
# https://docs.djangoproject.com/en/1.9/ref/settings/#databases

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': 'cepdb',
        'USER': 'cepdb',
        'PASSWORD': '599sOrjnP20[}a;f>JYHj/Znn-bn*V',
        'HOST': 'db-user.rc.fas.harvard.edu',
        'PORT': '5432',
    },

    'mongo': {
       'ENGINE': 'django.db.backends.dummy',
       'NAME': 'cep_syn',
       'USER': 'genAdmin',
       'PASSWORD': 'PceVb5LVJFjT',
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
