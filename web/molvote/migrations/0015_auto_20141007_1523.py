# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('molvote', '0014_auto_20141007_1519'),
    ]

    operations = [
        migrations.RenameField(
            model_name='vote',
            old_name='ballots',
            new_name='ballot',
        ),
    ]
