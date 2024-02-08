"""
Django settings for celeste project.

For more information on this file, see
https://docs.djangoproject.com/en/1.6/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/1.6/ref/settings/
"""

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
import os
import djcelery
import flower

from django.contrib import messages

djcelery.setup_loader()
BASE_DIR = os.path.dirname(os.path.dirname(__file__))

# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/1.6/howto/deployment/checklist/

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = '(*zdyo3k73l1ogo06irzvr$$orbivg#ijrv2$9zxa-e24rd0@h'

#CELERY SETUP
BROKER_URL = 'redis://localhost:6379/0'
CELERY_RESULT_BACKEND = 'redis://localhost:6379/0'

CELERYBEAT_SCHEDULER = "djcelery.schedulers.DatabaseScheduler"


# Application definition

INSTALLED_APPS = (
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'django.contrib.sites',
    'tastypie',
    'crispy_forms',
    'debug_toolbar',
    'django_extensions',
    'tastypie_mongoengine',
    'automate',
    'molvote',
    'utils',
    'djcelery',
    'datetimewidget',
    'bootstrap3',
    'registration',
    'invitekey',
    'dbbackup',
    'chemworkers',
    'rest_framework',
    'machinelearning'
)

MIDDLEWARE_CLASSES = (
    'debug_toolbar.middleware.DebugToolbarMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
    'django.middleware.locale.LocaleMiddleware'
)


ROOT_URLCONF = 'celeste.urls'

WSGI_APPLICATION = 'celeste.wsgi.application'



# Internationalization
# https://docs.djangoproject.com/en/1.6/topics/i18n/

LANGUAGE_CODE = 'en-us'

TIME_ZONE = 'America/New_York'

USE_I18N = True

USE_L10N = True

USE_TZ = True


# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/1.6/howto/static-files/

# set static to samsung/web/static
# this matches settings in nginx site settings file: docker-default
STATIC_ROOT = os.path.join(os.path.dirname(BASE_DIR), "static")
STATIC_URL = '/static/'

FLOWER_STATIC = os.path.join(os.path.dirname(flower.__file__),"static")

STATICFILES_DIRS = (
    os.path.join(os.path.dirname(BASE_DIR), "global_static"),
    FLOWER_STATIC
)


DEBUG_TOOLBAR_PATCH_SETTINGS = False
INTERNAL_IPS = ["127.0.0.1"]

CRISPY_TEMPLATE_PACK = 'bootstrap3'

### ALL AUTH

SITE_ID = 1

TEMPLATE_DIRS = (os.path.join(os.path.dirname(BASE_DIR), "templates"),)
TEMPLATE_CONTEXT_PROCESSORS = (
    "django.contrib.auth.context_processors.auth",
    "django.core.context_processors.debug",
    "django.core.context_processors.i18n",
    "django.core.context_processors.media",
    "django.core.context_processors.static",
    "django.core.context_processors.tz",
    "django.contrib.messages.context_processors.messages",

    # Required by allauth template tags
    "django.core.context_processors.request",

    # allauth specific context processors
    #"allauth.account.context_processors.account",
    #"allauth.socialaccount.context_processors.socialaccount",
)

AUTHENTICATION_BACKENDS = (
    # Needed to login by username in Django admin, regardless of `allauth`
    "django.contrib.auth.backends.ModelBackend",

)

LOGIN_REDIRECT_URL = "/"

MESSAGE_TAGS = {
    messages.ERROR: 'danger'
}

APPEND_SLASH = True


EMAIL_USE_TLS = True
EMAIL_HOST = 'smtp.mailgun.org'
EMAIL_HOST_USER = 'postmaster@sandbox69528288d63e430484ea0d9eaaeae380.mailgun.org'
EMAIL_HOST_PASSWORD = 'c3l3st3'
EMAIL_PORT = 587

ACCOUNT_ACTIVATION_DAYS=3

############### DB BACKUP SETTINGS #############
DBBACKUP_STORAGE = 'dbbackup.storage.dropbox_storage'
DBBACKUP_TOKENS_FILEPATH = os.path.join(BASE_DIR, "settings", "dropback_backup_tokens.txt")
DBBACKUP_DROPBOX_APP_KEY = '73oh4z5ls2z3ka6'
DBBACKUP_DROPBOX_APP_SECRET = 'uzsxzuzv3jlf290'

# DBBACKUP_STORAGE = 'dbbackup.storage.filesystem_storage'
# DBBACKUP_BACKUP_DIRECTORY = '/home/hirzel/db_backups/'

CACHES = {
    'default': {
        'BACKEND': 'redis_cache.RedisCache',
        'LOCATION': '127.0.0.1:6379',
        'OPTIONS': {
            'DB': 1
        }
    }

#     'default': {
#         'BACKEND': 'django.core.cache.backends.db.DatabaseCache',
#         'LOCATION': 'celeste_cache_table',
#     }
}


CHEM_WORKERS_LOG_PATH = "/var/log/chemworkers/"
CHEM_WORKERS_CONFIG_BASE = os.path.join(BASE_DIR, "..", "..", "workers", "configs")

IMG_CACHE_DIR = "/tmp"

REST_FRAMEWORK = {
    'DEFAULT_PERMISSION_CLASSES': ('rest_framework.permissions.IsAdminUser',),
    'PAGE_SIZE': 10,
    'DEFAULT_FILTER_BACKENDS': ('rest_framework.filters.DjangoFilterBackend',),
    'DEFAULT_PAGINATION_CLASS': 'rest_framework.pagination.LimitOffsetPagination'
}

NOTEBOOK_ARGUMENTS = [
    '--debug',
    '--config', '/home/blue/samsung/notebook/jupyter_notebook_config.py',
    '--notebook-dir', '/n/aagfs01/samsung/notebooks'
]

