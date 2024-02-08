# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('molvote', '0018_auto_20141205_1205'),
    ]

    operations = [
        migrations.AlterField(
            model_name='candidate',
            name='inchi_key',
            field=models.CharField(max_length=27),
        ),
    ]
