from .base import *


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

}

