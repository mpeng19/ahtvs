# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings
import datetime


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('molvote', '0008_auto_20141001_1133'),
    ]

    operations = [
        migrations.AddField(
            model_name='candidateset',
            name='announced',
            field=models.BooleanField(default=False),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='candidateset',
            name='creation_time',
            field=models.DateTimeField(null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='candidateset',
            name='creator',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL, null=True),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='ballot',
            name='creation_time',
            field=models.DateTimeField(default=datetime.datetime.now, null=True),
        ),
    ]
