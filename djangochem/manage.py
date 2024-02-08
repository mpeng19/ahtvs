#!/usr/bin/env python
import os
import sys
import logging

MY_DIR = os.path.dirname(__file__)
sys.path.append(os.path.abspath(os.path.join(MY_DIR, '.')))
sys.path.append(os.path.abspath(os.path.join(MY_DIR, '..')))

if __name__ == "__main__":
    if len(sys.argv) <= 1:
        print('No submmand provided.')
    if len(sys.argv) > 1 and sys.argv[1] == 'test':
        os.environ.setdefault("DJANGO_SETTINGS_MODULE", "djangochem.settings.test")
        logging.disable(logging.CRITICAL)
    else:
        os.environ.setdefault("DJANGO_SETTINGS_MODULE", "djangochem.settings.default")

    def handle_error(err, msg):
        if hasattr(err, 'msg'):
            err.msg += msg
        else:
            err.msg = msg
        raise err

    try:
        from django.core.management import execute_from_command_line
        from django.db.utils import OperationalError, ProgrammingError
    except ImportError as err:
        handle_error(err, "  Did you remember to run 'source activate a2g2'?")
    try:
        execute_from_command_line(sys.argv)
    except OperationalError as err:
        handle_error(err, "  Do you need to start your database server?")
    except ProgrammingError as err:
        handle_error(err, "  You may need to run 'django migrate'.")
