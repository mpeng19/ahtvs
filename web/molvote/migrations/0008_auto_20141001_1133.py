# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('molvote', '0007_auto_20141001_1108'),
    ]

    operations = [
        migrations.AddField(
            model_name='ballot',
            name='candidateset',
            field=models.ForeignKey(to='molvote.CandidateSet', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='ballot',
            name='name',
            field=models.CharField(default='noname', max_length=255),
            preserve_default=False,
        ),
        migrations.AlterField(
            model_name='candidate',
            name='calc_time',
            field=models.DateTimeField(null=True),
        ),
    ]
