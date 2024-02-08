# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('machinelearning', '0003_fingerprint'),
    ]

    operations = [
        migrations.AddField(
            model_name='neuralnet',
            name='fingerprint_name',
            field=models.CharField(max_length=32, null=True, blank=True),
        ),
        migrations.AlterField(
            model_name='fingerprint',
            name='name',
            field=models.CharField(max_length=32),
        ),
    ]
