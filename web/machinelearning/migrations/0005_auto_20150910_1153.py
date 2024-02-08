# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('molvote', '0027_auto_20150723_1229'),
        ('machinelearning', '0004_auto_20150910_1118'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='fingerprint',
            name='name',
        ),
        migrations.AddField(
            model_name='fingerprint',
            name='method',
            field=models.ForeignKey(to='molvote.Method', null=True),
        ),
        migrations.AlterField(
            model_name='fingerprint',
            name='candidate',
            field=models.ForeignKey(to='molvote.Candidate', null=True),
        ),
    ]
