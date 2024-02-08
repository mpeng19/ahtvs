# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('machinelearning', '0008_auto_20150911_0915'),
    ]

    operations = [
        migrations.AddField(
            model_name='neuralnet',
            name='batch_size',
            field=models.IntegerField(default=1),
        ),
        migrations.AddField(
            model_name='neuralnet',
            name='epochs',
            field=models.IntegerField(default=50),
        ),
    ]
