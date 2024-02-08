# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.utils.timezone


class Migration(migrations.Migration):

    dependencies = [
        ('molvote', '0027_auto_20150723_1229'),
        ('machinelearning', '0001_initial'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='neuralnet',
            name='candidateset',
        ),
        migrations.AddField(
            model_name='neuralnet',
            name='created_at',
            field=models.DateTimeField(default=django.utils.timezone.now),
        ),
        migrations.AddField(
            model_name='neuralnet',
            name='name',
            field=models.CharField(max_length=255, null=True, blank=True),
        ),
        migrations.AddField(
            model_name='neuralnet',
            name='properties',
            field=models.ManyToManyField(to='molvote.Property'),
        ),
    ]
