import sys

from base import *

# Database
# https://docs.djangoproject.com/en/1.6/ref/settings/#databases
from mongoengine import connect


MONGO_DBNAME = "test"
MONGO_HOST = "mongo"
MONGO_USER = ""
MONGO_PASSWORD = ""

if MONGO_USER or MONGO_PASSWORD:
    connect(MONGO_DBNAME, host=MONGO_HOST, username=MONGO_USER, password=MONGO_PASSWORD)
else:
    connect(MONGO_DBNAME, host=MONGO_HOST)

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.mysql',
        'NAME': 'test',  # Or path to database file if using sqlite3.
        # The following settings are not used with sqlite3:
        'USER': 'root',
        'PASSWORD': 'mysecretpassword',
        'HOST': 'mysql',
        'PORT': '3306',
        'OPTIONS': {
                'sql_mode': 'traditional',
            }
        },
    'mongo': {
       'ENGINE': 'django.db.backends.dummy',
       'NAME': 'test',
       'USER': '',
       'PASSWORD': '',
       'HOST': 'mongo',
       'PORT': '',
    }
}

if 'test' in sys.argv or 'test_coverage' in sys.argv:  # Covers regular testing and django-coverage
    DATABASES['default']['ENGINE'] = 'django.db.backends.sqlite3'


# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = True

TEMPLATE_DEBUG = True

ALLOWED_HOSTS = []


INSTALLED_APPS = INSTALLED_APPS + ('django_nose',)
# ####### LOGGING ################# #

LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
            'simple': {
                        'format': '%(levelname)s %(message)s'
                        },
        },
    'handlers': {
        'request': {
            'level': 'DEBUG',
            'class': 'logging.FileHandler',
            'filename': '/var/log/django/request.log',
        },
        'console': {
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


DUMP_CSV_PATH = os.path.join(BASE_DIR, "..", "..", "demos", "4.scatter_zoom")

DBBACKUP_SERVER_NAME = "celeste"
DBBACKUP_BACKUP_TEMP_DIRECTORY = "/tmp"

DBBACKUP_STORAGE = 'dbbackup.storage.filesystem_storage'
DBBACKUP_BACKUP_DIRECTORY = "/home/blue/database_files/"


TEST_RUNNER = 'django_nose.NoseTestSuiteRunner'


NOTEBOOK_ARGUMENTS = [
    '--debug',
    '--config', '/home/blue/samsung/notebook/jupyter_notebook_config.py',
    '--notebook-dir', '/home/blue/notebooks'
]

