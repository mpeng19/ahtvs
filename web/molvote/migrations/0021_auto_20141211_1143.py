# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('molvote', '0020_create_redoxpair'),
    ]

    operations = [
        migrations.CreateModel(
            name='RedoxPairOfPair',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('project', models.CharField(default=b'', max_length=255, blank=True)),
                ('released', models.BooleanField(default=False)),
                ('redox_potential', models.FloatField(null=True)),
                ('redox_potential_error', models.FloatField(null=True)),
                ('high_rp_pair', models.ForeignKey(related_name='high_rp_of_pair', to='molvote.RedoxPair')),
                ('low_rp_pair', models.ForeignKey(related_name='low_rp_of_pair', to='molvote.RedoxPair')),
            ],
            options={
                'abstract': False,
            },
            bases=(models.Model,),
        ),
        migrations.RemoveField(
            model_name='redoxpair',
            name='intermediate',
        ),
        migrations.AddField(
            model_name='redoxpair',
            name='is_minimum',
            field=models.BooleanField(default=False),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='redoxpair',
            name='hydrated_oxidized',
            field=models.ForeignKey(related_name='hyd_oxidized_of_redox_pair', to='molvote.Candidate', null=True),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='redoxpair',
            name='oxidized',
            field=models.ForeignKey(related_name='oxidized_of_redox_pair', to='molvote.Candidate'),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='redoxpair',
            name='reduced',
            field=models.ForeignKey(related_name='reduced_of_redox_pair', to='molvote.Candidate'),
            preserve_default=True,
        ),
    ]
