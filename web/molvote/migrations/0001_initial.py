# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Candidate',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('inchi_key', models.CharField(max_length=27)),
                ('smiles', models.CharField(max_length=1000)),
                ('absorption', models.FloatField()),
                ('splitting', models.FloatField()),
                ('strength', models.FloatField()),
                ('rate', models.FloatField()),
                ('homo', models.FloatField()),
                ('lumo', models.FloatField()),
                ('nicknames', models.CharField(max_length=1000)),
                ('date', models.DateField()),
                ('sascore', models.FloatField()),
                ('weight', models.FloatField()),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='CandidateSet',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=255)),
                ('members', models.ManyToManyField(to='molvote.Candidate')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.AddField(
            model_name='candidate',
            name='sets',
            field=models.ManyToManyField(to='molvote.CandidateSet'),
            preserve_default=True,
        ),
    ]
