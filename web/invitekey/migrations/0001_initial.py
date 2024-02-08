# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import invitekey.models


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Invitation',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('key', models.CharField(default=invitekey.models.make_key, unique=True, max_length=20)),
                ('project', models.CharField(max_length=50)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
    ]
