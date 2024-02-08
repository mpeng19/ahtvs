"""
WSGI config for celeste project.

It exposes the WSGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/1.6/howto/deployment/wsgi/
"""

import os
import socket
import sys

from web.celeste.settings.select_by_hostname import set_default_settings_with_hostname

set_default_settings_with_hostname()

sys.path.append(os.path.abspath(os.path.join(os.path.basename(__file__), "..", "..")))
sys.path.append("/opt/rdkit-Release_2015_03_1")

from django.core.wsgi import get_wsgi_application
application = get_wsgi_application()
