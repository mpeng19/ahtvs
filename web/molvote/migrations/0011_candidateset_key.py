# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('molvote', '0010_auto_20141002_1211'),
    ]

    operations = [
        migrations.AddField(
            model_name='candidateset',
            name='key',
            field=models.CharField(max_length=100, unique=True, null=True, blank=True),
            preserve_default=True,
        ),
    ]
