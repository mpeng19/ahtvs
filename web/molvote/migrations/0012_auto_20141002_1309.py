# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import molvote.models


class Migration(migrations.Migration):

    dependencies = [
        ('molvote', '0011_candidateset_key'),
    ]

    operations = [
        migrations.AlterField(
            model_name='candidateset',
            name='key',
            field=models.CharField(default=molvote.models.make_key, max_length=100, unique=True, null=True, blank=True),
        ),
    ]
