# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('molvote', '0017_auto_20141023_1822'),
    ]

    operations = [
        migrations.AddField(
            model_name='candidate',
            name='acceptor',
            field=models.ForeignKey(related_name=b'acceptor_for', to='molvote.Candidate', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='candidate',
            name='bridge1',
            field=models.ForeignKey(related_name=b'bridge1_for', to='molvote.Candidate', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='candidate',
            name='bridge2',
            field=models.ForeignKey(related_name=b'bridge2_for', to='molvote.Candidate', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='candidate',
            name='donor',
            field=models.ForeignKey(related_name=b'donor_for', to='molvote.Candidate', null=True),
            preserve_default=True,
        ),
    ]
