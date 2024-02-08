# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('molvote', '0012_auto_20141002_1309'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='candidateset',
            name='members',
        ),
    ]
