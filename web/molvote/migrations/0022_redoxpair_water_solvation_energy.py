# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('molvote', '0021_auto_20141211_1143'),
    ]

    operations = [
        migrations.AddField(
            model_name='redoxpair',
            name='water_solvation_energy',
            field=models.FloatField(null=True),
            preserve_default=True,
        ),
    ]
