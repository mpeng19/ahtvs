# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('molvote', '0016_ballot_close_time'),
    ]

    operations = [
        migrations.AlterField(
            model_name='ballot',
            name='close_time',
            field=models.DateTimeField(null=True, blank=True),
        ),
        migrations.AlterField(
            model_name='ballot',
            name='completion_time',
            field=models.DateTimeField(null=True, blank=True),
        ),
    ]
