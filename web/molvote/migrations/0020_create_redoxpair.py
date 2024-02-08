# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('molvote', '0019_auto_20141205_1243'),
    ]

    operations = [
        migrations.CreateModel(
            name='RedoxPair',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('project', models.CharField(default=b'', max_length=255, blank=True)),
                ('released', models.BooleanField(default=False)),
                ('log_hyd_constant', models.FloatField(null=True)),
                ('redox_potential', models.FloatField()),
                ('hydrated_oxidized', models.ForeignKey(related_name='hyd_oxidized_of_redox', to='molvote.Candidate')),
                ('intermediate', models.ForeignKey(related_name='intermediate_of_redox', to='molvote.Candidate')),
                ('oxidized', models.ForeignKey(related_name='oxidized_of_redox', to='molvote.Candidate')),
                ('reduced', models.ForeignKey(related_name='reduced_of_redox', to='molvote.Candidate')),
            ],
            options={
                'abstract': False,
            },
            bases=(models.Model,),
        ),
        migrations.AddField(
            model_name='candidate',
            name='water_solvation_energy',
            field=models.FloatField(null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='candidate',
            name='total_energy',
            field=models.FloatField(null=True),
            preserve_default=True,
        ),
         migrations.AddField(
            model_name='redoxpair',
            name='log_hyd_constant_error',
            field=models.FloatField(null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='redoxpair',
            name='redox_potential_error',
            field=models.FloatField(null=True),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='redoxpair',
            name='redox_potential',
            field=models.FloatField(null=True),
            preserve_default=True,
        ),
    ]
