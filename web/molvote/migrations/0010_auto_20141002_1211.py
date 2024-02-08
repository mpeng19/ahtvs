# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.utils.timezone


class Migration(migrations.Migration):

    dependencies = [
        ('molvote', '0009_auto_20141002_1209'),
    ]

    operations = [
        migrations.AlterField(
            model_name='ballot',
            name='creation_time',
            field=models.DateTimeField(default=django.utils.timezone.now, null=True),
        ),
    ]
