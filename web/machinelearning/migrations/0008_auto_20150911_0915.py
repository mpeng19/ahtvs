# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('machinelearning', '0007_neuralnet_training_size'),
    ]

    operations = [
        migrations.AddField(
            model_name='neuralnet',
            name='learn_rate',
            field=models.FloatField(default=0.001604),
        ),
        migrations.AddField(
            model_name='neuralnet',
            name='momentum',
            field=models.FloatField(default=0.98371),
        ),
    ]
