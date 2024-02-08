# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('chemworkers', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='worklocation',
            name='label',
            field=models.CharField(default='Test Location', max_length=255),
            preserve_default=False,
        ),
    ]
