# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('molvote', '0004_auto_20140919_1217'),
    ]

    operations = [
        migrations.AddField(
            model_name='candidate',
            name='project',
            field=models.CharField(default=b'', max_length=255, blank=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='candidate',
            name='released',
            field=models.BooleanField(default=False),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='candidateset',
            name='project',
            field=models.CharField(default=b'', max_length=255, blank=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='candidateset',
            name='released',
            field=models.BooleanField(default=False),
            preserve_default=True,
        ),
    ]
