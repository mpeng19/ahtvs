# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('molvote', '0027_auto_20150723_1229'),
    ]

    operations = [
        migrations.CreateModel(
            name='NeuralNet',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('weights', models.BinaryField()),
                ('mean', models.FloatField(null=True)),
                ('std', models.FloatField(null=True)),
                ('init_scale', models.FloatField(default=0.1)),
                ('num_inputs', models.IntegerField(default=512)),
                ('num_hidden', models.IntegerField(default=100)),
                ('prop_name', models.CharField(max_length=255)),
                ('candidateset', models.ForeignKey(to='molvote.CandidateSet')),
                ('method', models.ForeignKey(to='molvote.Method')),
            ],
        ),
    ]
