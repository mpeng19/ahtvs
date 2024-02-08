# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('molvote', '0027_auto_20150723_1229'),
        ('machinelearning', '0002_auto_20150908_1523'),
    ]

    operations = [
        migrations.CreateModel(
            name='Fingerprint',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('fingerprint', models.BinaryField()),
                ('name', models.CharField(max_length=255)),
                ('candidate', models.ForeignKey(to='molvote.Candidate')),
            ],
        ),
    ]
