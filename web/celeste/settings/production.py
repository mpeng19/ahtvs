from base import *

# Database
# https://docs.djangoproject.com/en/1.6/ref/settings/#databases
from mongoengine import connect

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = False
TEMPLATE_DEBUG = False

# this is only safe because we are running inside a protected environment where any other ports are not allowed
ALLOWED_HOSTS = ['*']

DUMP_CSV_PATH = os.path.join(BASE_DIR, "..","..","demos","4.scatter_zoom")

# DATABASES = {
#     'default': {
#         'ENGINE': 'django.db.backends.mysql',
#         'NAME': 'django_db',  # Or path to database file if using sqlite3.
#         # The following settings are not used with sqlite3:
#         'USER': 'blue',
#         'PASSWORD': 'c3l3st3',
#         'HOST': 'atlantic.chem.harvard.edu',
#         'PORT': '',
#     }
# }
DBBACKUP_SERVER_NAME = "celeste"
DBBACKUP_BACKUP_TEMP_DIRECTORY = "/n/aagfs01/samsung"

IMG_CACHE_DIR = "/home/hirzel/imgs"

DATABASES = {
   'default': {
       'ENGINE': 'django.db.backends.mysql',
       'NAME': 'celeste_db',
       'USER': 'celeste_db',
       'PASSWORD': 'sv5o:x}_KyB0TA6wP4&4',
       'HOST': 'db-user.rc.fas.harvard.edu',
       'PORT': '',
   },
   'mongo': {
       'ENGINE': 'django.db.backends.dummy',
       'NAME': 'samsung',
       'USER': 'blue',
       'PASSWORD': '2hXlaluy1ord',
       'HOST': 'molspace.rc.fas.harvard.edu',
       'PORT': '',
       'ALIAS': 'default',
       'PROJECTS': ['samsung']
   },
    'flow_batt': {
       'ENGINE': 'django.db.backends.dummy',
       'NAME': 'flow_batt',
       'USER': 'flow_batt',
       'PASSWORD': 'muchpotential!',
       'HOST': 'molspace.rc.fas.harvard.edu',
       'PORT': '',
       'ALIAS': 'flow_batt',
       'PROJECTS': ['flow_batt']
   }
}

mongo = DATABASES["mongo"]
connect(mongo["NAME"], alias=mongo["ALIAS"], host=mongo["HOST"], username=mongo["USER"], password=mongo["PASSWORD"])
fb = DATABASES["flow_batt"]
connect(fb["NAME"], alias=fb["ALIAS"], host=fb["HOST"], username=fb["USER"], password=fb["PASSWORD"])

######## LOGGING ##################

LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
            'simple': {
                'format': '%(levelname)s %(message)s'
            },
            'timestamped': {
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
            'formatter': 'timestamped'
        },
    },
    'loggers': {
        'django': {
            'handlers': ['error'],
            'propagate': True,
            'level': 'WARNING',
        },
        'django.request': {
            'handlers': ['request'],
            'level': 'WARNING',
            'propagate': True,
        },
    },
}

# this makes sure the DRF uses https with the HyperlinkedModelSerializer
SECURE_PROXY_SSL_HEADER = ('HTTP_X_FORWARDED_PROTO', 'https')
