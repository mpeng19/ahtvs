# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('molvote', '0024_candidate_cas'),
    ]

    operations = [
        migrations.AddField(
            model_name='redoxpair',
            name='michael_hydration_energy',
            field=models.FloatField(null=True),
            preserve_default=True,
        ),
    ]
