#!/usr/bin/env python

import sys

from celeste.settings.select_by_hostname import set_default_settings_with_hostname


set_default_settings_with_hostname()

from django.core.management import execute_from_command_line

execute_from_command_line(sys.argv)
