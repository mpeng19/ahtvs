from .base import *


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
       'NAME': 'test',
       'USER': '',
       'PASSWORD': '',
       'HOST': 'localhost',
       'PORT': '',
       'ALIAS': 'default',
       'PROJECTS': ['samsung']
    }
}

