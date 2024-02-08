# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.utils.timezone


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='ChemWorker',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(unique=True, max_length=255)),
                ('label', models.CharField(max_length=256)),
                ('config_path', models.CharField(max_length=512)),
                ('ignore_lock', models.BooleanField(default=False)),
                ('request_count_min', models.IntegerField(default=0)),
                ('request_count_max', models.IntegerField(default=10)),
                ('purge_mode', models.BooleanField(default=False)),
                ('min_priority', models.IntegerField(default=0)),
                ('max_priority', models.IntegerField(default=10)),
                ('project', models.CharField(max_length=256)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='WorkLocation',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('request_path', models.CharField(max_length=512)),
                ('complete_path', models.CharField(max_length=512)),
                ('archive_path', models.CharField(max_length=512)),
                ('error_path', models.CharField(max_length=512)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='WorkLog',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('timestamp', models.DateTimeField(default=django.utils.timezone.now)),
                ('jobs_requested', models.IntegerField(default=0)),
                ('jobs_completed', models.IntegerField(default=0)),
                ('jobs_failed', models.IntegerField(default=0)),
                ('chem_worker', models.ForeignKey(to='chemworkers.ChemWorker')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.AddField(
            model_name='chemworker',
            name='work_location',
            field=models.ForeignKey(to='chemworkers.WorkLocation'),
            preserve_default=True,
        ),
    ]
