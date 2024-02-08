# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('molvote', '0015_auto_20141007_1523'),
    ]

    operations = [
        migrations.AddField(
            model_name='ballot',
            name='close_time',
            field=models.DateTimeField(null=True),
            preserve_default=True,
        ),
    ]
