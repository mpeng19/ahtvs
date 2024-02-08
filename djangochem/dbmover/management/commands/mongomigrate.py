#!/usr/bin/env python

import pytz
import time
import logging

from rdkit.Chem import AllChem as Chem
import aag_python.molecular_storage.molecular_data_models as mdm
from django.contrib.auth.models import Group
from guardian.shortcuts import assign_perm
from pgmols.models import Mol, Calc, Geom, Method, Reaction
from jobs.models import Job, WorkBatch, JobConfig
from django.core.management.base import BaseCommand


PT = Chem.GetPeriodicTable()


PROP_MIGRATE = {"orbital_energy_list": "orbitals",
                "normal_mode_list": "normalmodes",
                "excited_states": "excitedstates",
                "velocity_list": "velocities"}


class Command(BaseCommand):
    help = 'migrate the connected mongo to connected postgres'

    def add_arguments(self, parser):
        parser.add_argument('input_project', type=str, default=None)
        parser.add_argument('output_project', type=str, default=None)
        parser.add_argument('-i', '--inchikey', type=str, default=None)
        parser.add_argument('-t', '--tag', type=str, default=None)
        parser.add_argument('-m', '--mol_only', action='store_true', default=False)
        parser.add_argument('-a', '--release_all', action='store_true', default=False)

    def handle(self, *args, **options):
        if options['input_project'] is None:
            raise Exception("a input_project must be specified")
        if options['output_project'] is None:
            raise Exception("a output_project must be specified")
        main(**options)


def build_props(mongo_calc):
    props = {}
    if mongo_calc.properties:
        for k, v in mongo_calc.properties._data.items():
            if v and not k.startswith('_'):
                new_key = k.replace('_', '')
                if new_key.endswith('list'):
                    new_key = new_key[:-4] + 's'
                props[new_key] = v
    for mongo_attr, pg_attr in PROP_MIGRATE.items():
        mongoval = getattr(mongo_calc, mongo_attr)
        if mongoval:
            props[pg_attr] = mongoval
    return props


def build_method(mongo_calc, job):
    if mongo_calc.theory is None:
        method, created = Method.objects.get_or_create(name='Unknown',
                                                       description='')
    else:
        description = mongo_calc.theory.theory_description
        if description.startswith("MMFF conformer"):
            description, params = description.split('.', 1)
            try:
                params = dict([line.split('=') for line in params.split(' ') if '=' in line])
                job.details['confgen_params'] = params
            except:
                job.details['confgen_params'] = params
            job.save()

        method, created = Method.objects.get_or_create(name=mongo_calc.theory.theory_level,
                                                       description=description)

    if created:
        logging.info("New method created {} ({})".format(method.name, method.description))
    return method


def build_geom(mongo_calc, method, job, parent_pgmol):
    geom = None
    details = {}
    if mongo_calc.coord_list:
        coords = []
        for atom_dict in mongo_calc.coord_list:
            row = [PT.GetAtomicNumber(str(atom_dict['element'])),
                   float(atom_dict['x']),
                   float(atom_dict['y']),
                   float(atom_dict['z'])]
            coords.append(row)
        if mongo_calc.meta_data.worker_name == 'conformer':
            for tag in mongo_calc.meta_data.tag_list:
                if "conf" in tag:
                    confnum = int(tag.split('_')[-1])
                    details['confnum'] = confnum
                    break
        geom = Geom.objects.create(method=method,
                                   parentjob=job,
                                   mol=parent_pgmol,
                                   xyz=coords,
                                   details=details)

    return geom


def ingest_calc(mongo_calc, parent_pgmol, group, parent_pggeom=None, workbatch=None):
    try:
        meta = mongo_calc.meta_data
    except AttributeError:
        try:
            mongo_calc = mdm.Calculation.objects.get(mongo_calc.id)
        except:
            return

    timeststamp = meta.document_creation_date.replace(tzinfo=pytz.UTC)
    details = dict(key=meta.key,
                   software=meta.program,
                   version=meta.version,
                   nodename=meta.host,
                   )
    # parent job
    if meta.worker_name is None:
        meta.worker_name = 'None'
    jobconfig, created = JobConfig.objects.get_or_create(name=meta.worker_name)
    if parent_pggeom:
        parent = parent_pggeom
    else:
        parent = parent_pgmol
    parent_job = Job.objects.create(status='done',
                                    duration=meta.wall_clock_time,
                                    cores=meta.num_procs,
                                    createtime=timeststamp,
                                    details=details,
                                    config=jobconfig,
                                    workbatch=workbatch,
                                    group=group,
                                    parent=parent)

    method = build_method(mongo_calc, parent_job)
    geom = build_geom(mongo_calc, method, parent_job, parent_pgmol)
    props = build_props(mongo_calc)

    pgcalc = Calc.objects.create(method=method,
                                 parentjob=parent_job,
                                 mol=parent_pgmol,
                                 props=props
                                 )

    # child jobs that have been ordered or are running now.
    if mongo_calc.meta_data.status.current_work_order:
        if mongo_calc.meta_data.status == 'clc_run':
            status = 'claimed'
        else:
            status = ''
        childjobconfig, created = JobConfig.objects.get_or_create(name=mongo_calc.meta_data.status.current_work_order)
        Job.objects.create(status=status,
                           config=childjobconfig,
                           group=group,
                           parent=geom)

    if geom:  # create a single point hanging off the new geom
        if parent_pggeom:
            geom.parents.add(parent_pggeom)
        pgcalc.geoms.add(geom)
    else:  # create just a single point
        if parent_pggeom:
            pgcalc.geoms.add(parent_pggeom)

    for mc in mongo_calc.child_calculation_list:
        ingest_calc(mc, parent_pgmol, group, geom, workbatch)


def ingest_mol(mongo_mol, group, workbatch, release_group, release_all, mol_only):
    if not hasattr(mongo_mol, "meta_data"):
        logging.error("corrupt molecule")
        return None
    pgmol, created = Mol.objects.get_or_create(smiles=mongo_mol.meta_data.smiles,
                                               group=group,
                                               inchikey=mongo_mol.meta_data.inchi_key)
    if created:
        pgmol.createtime = mongo_mol.meta_data.document_creation_date.replace(tzinfo=pytz.UTC)
        pgmol.tags = mongo_mol.meta_data.tag_list
        pgmol.nicknames = mongo_mol.nickname_list
        pgmol.mass = mongo_mol.mass
        pgmol.save()
        for linkage in mongo_mol.parent_list:
            if hasattr(linkage, "parent_molecule_list"):
                for mm in linkage.parent_molecule_list:
                    parent = ingest_mol(mm, group, workbatch, release_group, release_all, mol_only)
                    if parent is not None and parent != pgmol:
                        pgmol.parents.add(parent)

        if not mol_only:
            for mc in mongo_mol.calculation_list:
                ingest_calc(mc, pgmol, group, workbatch=workbatch)
        if mongo_mol.meta_data.status.current_work_order:
            if mongo_mol.meta_data.status == 'clc_run':
                status = 'claimed'
            else:
                status = ''
            config_name = mongo_mol.meta_data.status.current_work_order
            childjobconfig, created = JobConfig.objects.get_or_create(name=config_name)
            Job.objects.create(status=status,
                               config=childjobconfig,
                               group=group,
                               parent=pgmol)
    for child_mongomol in mongo_mol.child_list:
        for rtype in ['hydrate', 'bromide', 'michael_hyd', 'michael_hyd_R']:
            if rtype in child_mongomol.meta_data.tag_list:
                child_pgmol = ingest_mol(child_mongomol, group, workbatch, release_group, release_all, mol_only)
                if child_pgmol is not None:
                    reaction = Reaction(name=rtype)
                    reaction.save()
                    reaction.reactants.add(pgmol)
                    reaction.products.add(child_pgmol)
                    reaction.save()

    if mongo_mol.meta_data.status.released or release_all:
        assign_perm('view_mol', release_group, pgmol)

    return pgmol


def main(input_project, output_project, inchikey=None, release_all=False, mol_only=False, tag=None, **kwargs):
    starttime = time.time()
    release_group, created = Group.objects.get_or_create(name=output_project+"_nonadmin")
    workbatch, created = WorkBatch.objects.get_or_create(name='migration')
    if tag:
        query = mdm.Molecule.objects.filter(meta_data__project_name=input_project,
                                            meta_data__tag_list=tag).timeout(False).no_cache()
    elif inchikey:
        query = mdm.Molecule.objects.filter(meta_data__project_name=input_project,
                                            meta_data__inchi_key=inchikey).timeout(False).no_cache()
    else:
        query = mdm.Molecule.objects.filter(meta_data__project_name=input_project).timeout(False).no_cache()

    group, created = Group.objects.get_or_create(name=output_project)
    count = 0
    print(query.count(), 'molecules found')
    for mongo_mol in query:
        try:
            if ingest_mol(mongo_mol, group, workbatch, release_group, release_all, mol_only):
                count += 1
                if count % 100 == 0:
                    duration = time.time() - starttime
                    print("{} seconds.  {}/s".format(duration, count/duration))
        except:
            logging.exception(mongo_mol.meta_data.inchi_key)
            raise