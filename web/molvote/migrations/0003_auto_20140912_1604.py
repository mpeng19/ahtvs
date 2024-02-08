# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('molvote', '0002_auto_20140912_1031'),
    ]

    operations = [
        migrations.AlterField(
            model_name='candidate',
            name='date',
            field=models.DateTimeField(),
        ),
    ]
