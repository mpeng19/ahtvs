# -*- coding: utf-8 -*-
# Generated by Django 1.9.5 on 2016-04-12 17:43
from __future__ import unicode_literals

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('pgmols', '0003_auto_20160412_1742'),
    ]

    operations = [
        migrations.AlterIndexTogether(
            name='mol',
            index_together=set([]),
        ),
    ]
