"""
make new job objects in the DB
"""
import os
import importlib 
import json

from django.contrib.auth.models import User, Group

from pgmols.models import Geom, Mol, Calc
from jobs.models import Job, WorkBatch, JobConfig
from .dbinterface import postgresinterface 
from django.db import transaction

MAX_BEFORE_WARN = 10000

CLASS_FROM_STRING = {'Calc': Calc,
                     'Mol': Mol,
                     'Geom': Geom}

def create_workbatch(batch_name=None, username=None, comments=None):
    if username:
        user = User(name=username)
    else:
        user = None
    newbatch = WorkBatch(name=batch_name, user=user, comments=comments)
    newbatch.save()
    return newbatch

@transaction.atomic
def clear_error_status(project_name,
                       request_config,
                       limit=None):
    group = Group.objects.get(name=project_name)
    query = Job.objects.select_for_update().filter(group=group, status=postgresinterface.ERROR, config__name=request_config)
    if limit is not None:
        query = query.limit(limit)
    return query.update(status=postgresinterface.CLEAR)

#@transaction.atomic
def request_jobs(project,
                 requested_config_name,
                 parent_config_name=None,
                 limit=None,
                 workbatch=None,
                 details={},
                 avoid_dupes=True,
                 tags=[],
                 force=False,
                 request_spec=None,
                 dry_run=None):
    config = JobConfig.objects.get(name=requested_config_name)
    group = Group.objects.get(name=project)
    request_builder_path = os.path.join(config.dir, 'job_request_builder.py')

    with transaction.atomic():
        if os.path.exists(request_builder_path):
            job_dicts = build_job_dicts_from_config_builder(
                config=config,
                builder_path=request_builder_path,
                request_spec=request_spec)
        else:
            job_dicts = build_job_dicts_on_objects(
                config=config,
                group=group,
                tags=tags,
                limit=limit,
                force=force,
                parent_config_name=parent_config_name,
                workbatch=workbatch,
                avoid_dupes=avoid_dupes)

        if dry_run:
            return job_dicts
        else:
            return create_jobs_from_dicts(
                job_dicts=job_dicts,
                workbatch=workbatch,
                config=config,
                group=group,
                detail_overrides=details)


def build_job_dicts_from_config_builder(config=None, builder_path=None,
                                        request_spec=None):
    module_name = config.name + '.' + 'job_request_builder'
    spec = importlib.util.spec_from_file_location(module_name, builder_path)
    builder = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(builder)
    return builder.build_job_dicts(request_spec=request_spec)


def build_job_dicts_on_objects(config=None,
                               group=None,
                               tags=None,
                               parent_config_name=None,
                               limit=None,
                               force=None,
                               workbatch=None,
                               avoid_dupes=None):
    parent_model = CLASS_FROM_STRING[config.parent_class_name]
    if parent_model == Mol:
        parent_query = Mol.objects.filter(group=group)
        if tags:
            parent_query = parent_query.filter(tags__contains=tags)
    else:
        if parent_config_name is None and not force:
            raise Exception("Please specify a parent_config for this job"
                            " request. (Use force to bypass this warning"
                            " if you really know what you are doing)")
        parent_query = parent_model.objects.filter(mol__group=group)
        if tags:
            parent_query = parent_query.filter(mol__tags__contains=tags)
    if parent_config_name is not None:
        parent_query = parent_query.filter(
            parentjob__config__name=parent_config_name)
    if workbatch:
        parent_query = parent_query.filter(parentjob__workbatch=workbatch)
    if avoid_dupes:  # prevent ordering the same job twice
        parent_query = parent_query.exclude(
            childjobs__config__name=config.name)
    if parent_query.count() > MAX_BEFORE_WARN and not force and limit is None:
        raise Exception("Sorry, you are attempt to create too many ({}) new"
                        " jobs. Use force to bypass this warning.".format(
                            parent_query.count()))
    if limit:
        parent_query = parent_query[:limit]
    job_dicts = []
    for parent in parent_query:
        job_dict = {'parent': parent}
        if parent.details:
            job_dict['details'] = {**parent.details}
        job_dicts.append(job_dict)
    return job_dicts


def create_jobs_from_dicts(job_dicts=None, workbatch=None, config=None,
                           group=None, detail_overrides=None):
    config_details = get_config_details(config) or {}
    detail_overrides = detail_overrides or {}
    jobs = []
    for job_dict in job_dicts:
        job_dict['details'] = {**job_dict.get('details', {}), **config_details,
                               **detail_overrides}
        jobs.append(Job(
            workbatch=workbatch,
            config=config,
            group=group,
            **job_dict))
    return Job.objects.bulk_create(jobs)

def get_config_details(config):
    if config.configpath:
        config_dict = json.load(open(config.configpath))
        if "default_details_filename" in config_dict:
            defaults_filename = config_dict["default_details_filename"]
            defaults_path = os.path.join(os.path.dirname(config.configpath), defaults_filename)
            return json.load(open(defaults_path))


# LEGACY, still used by at least wcgrid tests :(
def request_job_on_objects(parent_query,
                           requested_config,
                           group_name,
                           limit=None,
                           workbatch=None,
                           details={},
                           inherit_details=True):
    if limit:
        parent_query = parent_query[:limit]

    config = JobConfig.objects.get(name=requested_config)
    group = Group.objects.get(name=group_name)
    common = dict(workbatch=workbatch, config=config, group=group)

    defaults_dict = {}
    if config.configpath:
        config_dict = json.load(open(config.configpath))
        if "default_details_filename" in config_dict:
            defaults_filename = config_dict["default_details_filename"]
            defaults_path = os.path.join(os.path.dirname(config.configpath), defaults_filename)
            defaults_dict = json.load(open(defaults_path))

    def join_details(parent):
        if inherit_details and parent.details is not None:
            # parent details will be overwritten by details with same key
            return {**parent.details, **defaults_dict, **details}
        else:
            return {**defaults_dict, **details}

    new_jobs = [Job(parent=parent, details=join_details(parent), **common) for parent in parent_query]
    return Job.objects.bulk_create(new_jobs)
