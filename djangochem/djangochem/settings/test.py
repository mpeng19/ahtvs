import logging
from .base import *


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
       'NAME': 'unittest',
       'USER': '',
       'PASSWORD': '',
       'HOST': 'localhost',
       'PORT': '',
       'ALIAS': 'default',
       'PROJECTS': ['samsung']
    }
}

try:
    import mongoengine
    try:
        mongo = DATABASES["mongo"]
        mongoengine.connect(mongo["NAME"],
                            alias=mongo['ALIAS'],
                            host=mongo["HOST"],
                            username=mongo["USER"],
                            password=mongo["PASSWORD"])
    except mongoengine.connection.ConnectionError:
        logging.warning("No mongo database")
except ImportError:
    logging.warning("No Mongoengine")


LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'handlers': {
        'console': {
            'level': 'WARNING',
            'class': 'logging.StreamHandler',
        },
    },
    'loggers': {
        'django': {
            'handlers': ['console'],
            'level': 'DEBUG',
            'propagate': True,
        },
    },
}
