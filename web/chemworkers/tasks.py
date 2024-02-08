from __future__ import absolute_import
import logging
import os

from celery import task
from celery.utils.log import get_task_logger
from django.core.cache import cache
from django.conf import settings

from workers.worker import Worker, setup_log

from .models import ChemWorker, WorkLog

logger = logging.getLogger(__name__)

from aag_python.molecular_storage import molecular_data_models as mdm
from mongoengine.context_managers import switch_db

LOG_PATH = getattr(settings, "CHEM_WORKERS_LOG_PATH", "/var/log")
CONFIG_BASE = getattr(settings, "CHEM_WORKERS_CONFIG_BASE")

LOCK_EXPIRE = 60 * 15 # Lock expires in 15 minutes

## Use django database config to associate different projects with mongo databases
DB_ALIAS_LOOKUP = {}
for db_config in settings.DATABASES.values():
    for project_name in db_config.get("PROJECTS", []):
        DB_ALIAS_LOOKUP[project_name] = db_config["ALIAS"]


@task(name="chemworkers.do_one_round", ignore_result=True,
      store_errors_even_if_ignored=True)
def do_one_round(worker_name):
    logger = setup_log(worker_name, os.path.join(LOG_PATH, worker_name + ".log"))
    lock_id = 'lock-{}'.format(worker_name)
    # cache.add fails if if the key already exists

    def acquire_lock():
        return cache.add(lock_id, 'true', LOCK_EXPIRE)
    # memcache delete is very slow, but we have to use it to take
    # advantage of using add() for atomic locking

    def release_lock():
        cache.delete(lock_id)

    if acquire_lock():
        try:
            try:
                chem_worker = ChemWorker.objects.get(name=worker_name)
            except:
                logger.error("Could not find ChemWorker {}".format(worker_name))
                raise
            config_path = os.path.join(CONFIG_BASE, chem_worker.config_path)
            w = Worker(config_path=config_path,
                       inbox_job_dir=chem_worker.work_location.request_path,
                       completed_job_dir=chem_worker.work_location.complete_path,
                       archive_dir=chem_worker.work_location.archive_path,
                       error_dir=chem_worker.work_location.error_path,
                       ignore_lock=chem_worker.ignore_lock,
                       inbox_dir_count_limit=chem_worker.request_count_max,
                       inbox_trigger_count=chem_worker.request_count_min,
                       batch_count=chem_worker.batch_size,
                       purge_mode=chem_worker.purge_mode,
                       min_priority=chem_worker.min_priority,
                       project=chem_worker.project,
                       logger=logger
                       )
            with switch_db(mdm.Calculation, DB_ALIAS_LOOKUP[chem_worker.project]):
                with switch_db(mdm.Molecule, DB_ALIAS_LOOKUP[chem_worker.project]):
                    created_count, completed_count, error_count = w.run_once()

            if any([created_count, completed_count, error_count]):
                log = WorkLog(chem_worker=chem_worker,
                              jobs_requested=created_count,
                              jobs_completed=completed_count,
                              jobs_failed=error_count)
                log.save()
            else:
                logger.info("None created or completed.  No errors.")
        finally:
            release_lock()
    else:
        logger.warning("last task [{}] unfinished.  Skipping work round.".format(worker_name))