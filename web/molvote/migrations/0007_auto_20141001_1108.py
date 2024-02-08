# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('molvote', '0006_auto_20141001_0934'),
    ]

    operations = [
        migrations.CreateModel(
            name='Ballot',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('creation_time', models.DateTimeField(null=True)),
                ('completion_time', models.DateTimeField(null=True)),
                ('voter', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Comment',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('time', models.DateTimeField()),
                ('text', models.CharField(max_length=3000)),
                ('candidate', models.ForeignKey(to='molvote.Candidate')),
                ('user', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Vote',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('time', models.DateTimeField(null=True)),
                ('rating', models.IntegerField(null=True)),
                ('candidate', models.ForeignKey(to='molvote.Candidate')),
                ('voter', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.AddField(
            model_name='ballot',
            name='votes',
            field=models.ManyToManyField(to='molvote.Vote'),
            preserve_default=True,
        ),
        migrations.RenameField(
            model_name='candidate',
            old_name='date',
            new_name='calc_time'
        )
    ]
