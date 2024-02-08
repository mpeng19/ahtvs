# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('chemworkers', '0002_worklocation_label'),
    ]

    operations = [
        migrations.AddField(
            model_name='chemworker',
            name='batch_size',
            field=models.IntegerField(default=1),
            preserve_default=True,
        ),
    ]
