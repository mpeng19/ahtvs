# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('molvote', '0013_remove_candidateset_members'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='ballot',
            name='votes',
        ),
        migrations.RemoveField(
            model_name='vote',
            name='voter',
        ),
        migrations.AddField(
            model_name='vote',
            name='ballots',
            field=models.ForeignKey(to='molvote.Ballot', null=True),
            preserve_default=True,
        ),
    ]
