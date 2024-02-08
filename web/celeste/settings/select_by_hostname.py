import os
import socket 


def set_default_settings_with_hostname():
	hostname = socket.gethostname()
	
	if hostname == "atlantic":
	    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "web.celeste.settings.atlantic")
	elif hostname == "test":
	    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "web.celeste.settings.dev")
	else:
	    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "web.celeste.settings.production")