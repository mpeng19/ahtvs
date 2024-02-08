# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('molvote', '0003_auto_20140912_1604'),
    ]

    operations = [
        migrations.AlterField(
            model_name='candidate',
            name='absorption',
            field=models.FloatField(null=True),
        ),
        migrations.AlterField(
            model_name='candidate',
            name='homo',
            field=models.FloatField(null=True),
        ),
        migrations.AlterField(
            model_name='candidate',
            name='lumo',
            field=models.FloatField(null=True),
        ),
        migrations.AlterField(
            model_name='candidate',
            name='rate',
            field=models.FloatField(null=True),
        ),
        migrations.AlterField(
            model_name='candidate',
            name='sascore',
            field=models.FloatField(null=True),
        ),
        migrations.AlterField(
            model_name='candidate',
            name='splitting',
            field=models.FloatField(null=True),
        ),
        migrations.AlterField(
            model_name='candidate',
            name='strength',
            field=models.FloatField(null=True),
        ),
        migrations.AlterField(
            model_name='candidate',
            name='weight',
            field=models.FloatField(null=True),
        ),
    ]
