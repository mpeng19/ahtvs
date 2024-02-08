# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('machinelearning', '0009_auto_20150911_0928'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='neuralnet',
            name='properties',
        ),
        migrations.AlterField(
            model_name='neuralnet',
            name='batch_size',
            field=models.IntegerField(default=256),
        ),
    ]
