# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('molvote', '0022_redoxpair_water_solvation_energy'),
    ]

    operations = [
        migrations.AddField(
            model_name='redoxpair',
            name='substituent_count',
            field=models.IntegerField(null=True),
            preserve_default=True,
        ),
    ]
