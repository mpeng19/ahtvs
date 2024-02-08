from base import *

# Database
# https://docs.djangoproject.com/en/1.6/ref/settings/#databases
from mongoengine import connect

MONGO_DBNAME = "test"
MONGO_HOST = "localhost"

DBBACKUP_SERVER_NAME = "atlantic"

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.mysql',
        'NAME': 'django_db',
        'USER': 'blue',
        'PASSWORD': 'c3l3st3',
        'HOST': 'atlantic.chem.harvard.edu',
        'PORT': '',
    },
    'mongo': {
        'ENGINE': 'django.db.backends.dummy',
        'NAME': 'test',  # db name
        'USER': '',
        'PASSWORD': '',
        'HOST': 'localhost',
        'PORT': '',
    }
}

mongo = DATABASES["mongo"]
connect(mongo["NAME"], host=mongo["HOST"], username=mongo["USER"], password=mongo["PASSWORD"])

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = True

TEMPLATE_DEBUG = True

ALLOWED_HOSTS = []



######## LOGGING ##################

LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
            'simple': {
            'format': '%(asctime)s %(levelname)s %(message)s'
        },
        },
    'handlers': {
        'request': {
            'level': 'DEBUG',
            'class': 'logging.FileHandler',
            'filename': '/var/log/django/request.log',
        },
        'console':{
            'level': 'DEBUG',
            'class': 'logging.StreamHandler',
            'formatter': 'simple'
        },
        'error': {
            'level': 'DEBUG',
            'class': 'logging.FileHandler',
            'filename': '/var/log/django/error.log',
        },
    },
    'loggers': {
        'django': {
            'handlers': ['console', 'error'],
            'propagate': True,
            'level': 'INFO',
        },
        'django.request': {
            'handlers': ['request'],
            'level': 'DEBUG',
            'propagate': True,
        },
    },
}


DUMP_CSV_PATH = os.path.join(BASE_DIR, "..","..","demos","4.scatter_zoom")
