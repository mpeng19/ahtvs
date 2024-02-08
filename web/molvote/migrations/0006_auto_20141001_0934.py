# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('molvote', '0005_auto_20140925_1108'),
    ]

    operations = [
        migrations.AlterField(
            model_name='candidate',
            name='nicknames',
            field=models.CharField(default=b'', max_length=1000, blank=True),
        ),
    ]
