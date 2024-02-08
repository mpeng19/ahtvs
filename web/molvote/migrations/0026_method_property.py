# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('molvote', '0025_redoxpair_michael_hydration_energy'),
    ]

    operations = [
        migrations.CreateModel(
            name='Method',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('geom_theory', models.CharField(max_length=12, choices=[(b'mm', b'Molecular Mechanics'), (b'se', b'Semiempirical'), (b'dft', b'DFT'), (b'wf', b'Wave Function')])),
                ('geom_level', models.CharField(max_length=12, choices=[(b'DFT levels', [(b'gs', b'Ground State'), (b'rpa', b'RPA'), (b'tda', b'TDA')]), (b'MM Forcefields', [(b'mmff94', b'MMFF94'), (b'uff', b'UFF')]), (b'Semiempirical method', [(b'pm7', b'PM7')])])),
                ('geom_basis', models.CharField(max_length=12, choices=[(b'631gs', b'6-31G*'), (b'6311pgdp', b'6-311+Gdp')])),
                ('geom_functional', models.CharField(max_length=12, choices=[(b'b3lyp', b'B3LYP'), (b'm062x', b'M062X'), (b'wb97xd', b'Omega-B97xD'), (b'pbe', b'PBE')])),
                ('geom_solvation', models.CharField(max_length=12, choices=[(b'iefpcm', b'IEFPCM'), (b'cpcm', b'CPCM'), (b'cosmo', b'COSMO')])),
                ('geom_solvent', models.CharField(max_length=12, choices=[(b'water', b'Water'), (b'toluene', b'Toluene')])),
                ('geom_shell', models.CharField(max_length=12, choices=[(b'r', b'Restricted'), (b'u', b'Unrestricted'), (b'ro', b'Openshell')])),
                ('geom_multiplicity', models.IntegerField(choices=[(1, b'singlet'), (2, b'doublet'), (3, b'triplet')])),
                ('sp_theory', models.CharField(max_length=12, choices=[(b'mm', b'Molecular Mechanics'), (b'se', b'Semiempirical'), (b'dft', b'DFT'), (b'wf', b'Wave Function')])),
                ('sp_level', models.CharField(max_length=12, choices=[(b'DFT levels', [(b'gs', b'Ground State'), (b'rpa', b'RPA'), (b'tda', b'TDA')]), (b'MM Forcefields', [(b'mmff94', b'MMFF94'), (b'uff', b'UFF')]), (b'Semiempirical method', [(b'pm7', b'PM7')])])),
                ('sp_basis', models.CharField(max_length=12, choices=[(b'631gs', b'6-31G*'), (b'6311pgdp', b'6-311+Gdp')])),
                ('sp_functional', models.CharField(max_length=12, choices=[(b'b3lyp', b'B3LYP'), (b'm062x', b'M062X'), (b'wb97xd', b'Omega-B97xD'), (b'pbe', b'PBE')])),
                ('sp_solvation', models.CharField(max_length=12, choices=[(b'iefpcm', b'IEFPCM'), (b'cpcm', b'CPCM'), (b'cosmo', b'COSMO')])),
                ('sp_solvent', models.CharField(max_length=12, choices=[(b'water', b'Water'), (b'toluene', b'Toluene')])),
                ('sp_shell', models.CharField(max_length=12, choices=[(b'r', b'Restricted'), (b'u', b'Unrestricted'), (b'ro', b'Openshell')])),
                ('sp_multiplicity', models.IntegerField(choices=[(1, b'singlet'), (2, b'doublet'), (3, b'triplet')])),
                ('name', models.CharField(unique=True, max_length=255)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Property',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=255)),
                ('value', models.FloatField()),
                ('candidate', models.ForeignKey(to='molvote.Candidate')),
                ('method', models.ForeignKey(to='molvote.Method')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
    ]
