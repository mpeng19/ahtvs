# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('machinelearning', '0010_auto_20150911_1655'),
    ]

    operations = [
        migrations.AddField(
            model_name='neuralnet',
            name='l2_reg',
            field=models.FloatField(default=0.0),
        ),
    ]
