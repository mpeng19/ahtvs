# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('invitekey', '0002_auto_20141017_1545'),
    ]

    operations = [
        migrations.RenameField(
            model_name='invitation',
            old_name='project',
            new_name='group',
        ),
    ]
