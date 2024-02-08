# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('molvote', '0026_method_property'),
    ]

    operations = [
        migrations.AlterField(
            model_name='method',
            name='geom_basis',
            field=models.CharField(blank=True, max_length=12, choices=[(b'631gs', b'6-31G*'), (b'6311pgdp', b'6-311+Gdp')]),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='method',
            name='geom_functional',
            field=models.CharField(blank=True, max_length=12, choices=[(b'b3lyp', b'B3LYP'), (b'm062x', b'M062X'), (b'wb97xd', b'Omega-B97xD'), (b'pbe', b'PBE')]),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='method',
            name='geom_level',
            field=models.CharField(blank=True, max_length=12, choices=[(b'DFT levels', [(b'gs', b'Ground State'), (b'rpa', b'RPA'), (b'tda', b'TDA')]), (b'MM Forcefields', [(b'mmff94', b'MMFF94'), (b'uff', b'UFF')]), (b'Semiempirical method', [(b'pm7', b'PM7')])]),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='method',
            name='geom_multiplicity',
            field=models.IntegerField(null=True, choices=[(1, b'singlet'), (2, b'doublet'), (3, b'triplet')]),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='method',
            name='geom_shell',
            field=models.CharField(blank=True, max_length=12, choices=[(b'r', b'Restricted'), (b'u', b'Unrestricted'), (b'ro', b'Openshell')]),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='method',
            name='geom_solvation',
            field=models.CharField(blank=True, max_length=12, choices=[(b'iefpcm', b'IEFPCM'), (b'cpcm', b'CPCM'), (b'cosmo', b'COSMO')]),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='method',
            name='geom_solvent',
            field=models.CharField(blank=True, max_length=12, choices=[(b'water', b'Water'), (b'toluene', b'Toluene')]),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='method',
            name='geom_theory',
            field=models.CharField(blank=True, max_length=12, choices=[(b'mm', b'Molecular Mechanics'), (b'se', b'Semiempirical'), (b'dft', b'DFT'), (b'wf', b'Wave Function')]),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='method',
            name='sp_basis',
            field=models.CharField(blank=True, max_length=12, choices=[(b'631gs', b'6-31G*'), (b'6311pgdp', b'6-311+Gdp')]),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='method',
            name='sp_functional',
            field=models.CharField(blank=True, max_length=12, choices=[(b'b3lyp', b'B3LYP'), (b'm062x', b'M062X'), (b'wb97xd', b'Omega-B97xD'), (b'pbe', b'PBE')]),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='method',
            name='sp_level',
            field=models.CharField(blank=True, max_length=12, choices=[(b'DFT levels', [(b'gs', b'Ground State'), (b'rpa', b'RPA'), (b'tda', b'TDA')]), (b'MM Forcefields', [(b'mmff94', b'MMFF94'), (b'uff', b'UFF')]), (b'Semiempirical method', [(b'pm7', b'PM7')])]),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='method',
            name='sp_multiplicity',
            field=models.IntegerField(null=True, choices=[(1, b'singlet'), (2, b'doublet'), (3, b'triplet')]),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='method',
            name='sp_shell',
            field=models.CharField(blank=True, max_length=12, choices=[(b'r', b'Restricted'), (b'u', b'Unrestricted'), (b'ro', b'Openshell')]),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='method',
            name='sp_solvation',
            field=models.CharField(blank=True, max_length=12, choices=[(b'iefpcm', b'IEFPCM'), (b'cpcm', b'CPCM'), (b'cosmo', b'COSMO')]),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='method',
            name='sp_solvent',
            field=models.CharField(blank=True, max_length=12, choices=[(b'water', b'Water'), (b'toluene', b'Toluene')]),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='method',
            name='sp_theory',
            field=models.CharField(blank=True, max_length=12, choices=[(b'mm', b'Molecular Mechanics'), (b'se', b'Semiempirical'), (b'dft', b'DFT'), (b'wf', b'Wave Function')]),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='property',
            name='candidate',
            field=models.ForeignKey(related_name='properties', to='molvote.Candidate'),
            preserve_default=True,
        ),
    ]
