# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('molvote', '0023_redoxpair_substituent_count'),
    ]

    operations = [
        migrations.AddField(
            model_name='candidate',
            name='cas',
            field=models.CharField(max_length=12, null=True, blank=True),
            preserve_default=True,
        ),
    ]
