# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('molvote', '0027_auto_20150723_1229'),
        ('machinelearning', '0005_auto_20150910_1153'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='neuralnet',
            name='fingerprint_name',
        ),
        migrations.RemoveField(
            model_name='neuralnet',
            name='method',
        ),
        migrations.AddField(
            model_name='neuralnet',
            name='fp_method',
            field=models.ForeignKey(related_name='fp_method', to='molvote.Method', null=True),
        ),
        migrations.AddField(
            model_name='neuralnet',
            name='qc_method',
            field=models.ForeignKey(related_name='qc_method', to='molvote.Method', null=True),
        ),
        migrations.AlterField(
            model_name='neuralnet',
            name='prop_name',
            field=models.CharField(max_length=255, null=True, blank=True),
        ),
    ]
