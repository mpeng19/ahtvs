# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('invitekey', '0003_auto_20141017_1717'),
    ]

    operations = [
        migrations.AddField(
            model_name='invitation',
            name='email',
            field=models.CharField(default=b'', max_length=100),
            preserve_default=True,
        ),
    ]
