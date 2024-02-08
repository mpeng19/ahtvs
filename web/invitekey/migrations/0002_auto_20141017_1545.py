# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import invitekey.models


class Migration(migrations.Migration):

    dependencies = [
        ('invitekey', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='invitation',
            name='key',
            field=models.CharField(default=invitekey.models.make_key, unique=True, max_length=50),
        ),
    ]
