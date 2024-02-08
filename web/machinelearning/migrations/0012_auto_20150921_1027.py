# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('machinelearning', '0011_neuralnet_l2_reg'),
    ]

    operations = [
        migrations.RenameField(
            model_name='neuralnet',
            old_name='num_hidden',
            new_name='hidden_size',
        ),
        migrations.RenameField(
            model_name='neuralnet',
            old_name='num_inputs',
            new_name='input_size',
        ),
        migrations.AddField(
            model_name='neuralnet',
            name='hidden_layers',
            field=models.IntegerField(default=1),
        ),
    ]
