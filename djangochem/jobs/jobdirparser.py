import logging
import tempfile
import traceback
import importlib.util
import os
import json
import sys

import importlib

from .storage.filestorage import ReadDir
from .dbinterface import dbinterface as dbi

DEFAULT_PRIORITY = 5
DEFAULT_ERROR_FILE_NAME = "job_error.txt"
DEFAULT_JOB_INFO_FILENAME = "job_info.json"
DEFAULT_JOB_SCRIPT_FILENAME = "job.sh"
DEFAULT_TEMPLATE_PATH = "templates"
DEFAULT_PERIOD = 10
DEFAULT_LOG = "worker_error.log"
INBOX_DIR_COUNT_LIMIT = 10
INBOX_DIR_COUNT_TRIGGER = 0

DEFAULT_MIN_PRIORITY = 0

# PRIORITY_RANGE is the max possible range of priorities, tied to the priority system used in the database.
PRIORITY_RANGE = (0, 9)
TEMP_FILE_PREFIX = "temp_"

EXCEPTION_FILENAME = "job_loader_exception.txt"

logger = logging.getLogger('jobdirparser')

class BatchJobFound(UserWarning):
    pass


class ParseError(Exception):
    pass


class JobDirParser(object):

    """ Creates database objects from completed job dirs"""

    def __init__(self,
                 name,
                 project,
                 loader_path,
                 loader_module_name,
                 db_interface,
                 force=False,
                 storage_class=ReadDir,
                 storage_kwargs={}):
        self.name = name
        self.db_interface = db_interface
        self.project = project
        self.storage_class = storage_class
        self.storage_kwargs = storage_kwargs
        self.force = force
        try:
            module_path = os.path.join(loader_path, loader_module_name + ".py")
            spec = importlib.util.spec_from_file_location(loader_module_name,
                                                          module_path)
            self.loader_module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(self.loader_module)
        except Exception as e:
            print("Unable to load parser module for '%s': %s" % (
                name, e), file=sys.stderr)

    def update_db_for_completed_jobs(self, completed_job_dir, limit=None):
        """
        wrapper around iter_completed_jobs

        works as a single call instead of generator and logs results
        """
        for result, error in self.iter_completed_jobs(completed_job_dir, limit):
            if error:
                logger.warning(str(error))
            else:
                logger.info("parsed " + result)

    def iter_completed_jobs(self, completed_job_dir, limit=None,
                            lockfile=None):
        """
        generator to iterate through parsing of completed jobs.

        return tuple with the parsed job path and any exceptions that occurered
        """
        storage = self.storage_class(completed_job_dir, self.name, **self.storage_kwargs)
        for jd in storage.list_dir():
            if storage.dir_exists(jd):
                if lockfile and self.dir_has_lockfile(
                    storage=storage, dir=jd, lockfile=lockfile): continue
                batchjob = False
                try:
                    job_path = self.update_db_for_job(jd)
                    yield (job_path, None)
                except BatchJobFound:
                    batchjob = True
                except Warning as err:
                    yield (None, err)
                if batchjob:
                    batchdir = jd
                    storage = self.storage_class(batchdir, self.name, **self.storage_kwargs)
                    for subdir in storage.list_dir(dir_only=True):
                        try:
                            yield (self.update_db_for_job(subdir), None)
                        except Warning as err:
                            yield (None, err)
                    # clean up batch directory
                    batchdir = self.storage_class(jd, self.name, **self.storage_kwargs)
                    batchdir.archive()

            else:
                yield (None, ParseError("parser asked to handle file as job dir: {}".format(jd)))

    def dir_has_lockfile(self, storage=None, dir=None, lockfile=None):
        lockfile_path = os.path.join(dir, lockfile)
        return storage.file_exists(lockfile_path)

    def _verify_job(self, info, storage):
        # verify the jobdir config name matches the config we are working on
        jd_config_name = info.get(dbi.WORKER_NAME_KEY)
        if jd_config_name != self.name:
            msg = "Parser config {} does not match job dir creator {}"
            raise UserWarning(msg.format(self.name, jd_config_name))
        # verify this jobdir is for the project we are working on
        jd_project_name = info.get(dbi.PROJECT_KEY)
        if jd_project_name != self.project:
            msg = "Parser project {} does not match job dir project {}"
            raise UserWarning(msg.format(self.project, jd_project_name))
        # verify it's not a batch job
        if info.get(dbi.BATCH_KEY, False):
            raise BatchJobFound()
        # verify there is no error file
        if storage.file_exists(DEFAULT_ERROR_FILE_NAME):
            raise Exception("Error file exists for job dir {}".format(storage))
        # make sure an job exists in the db before bothering to parse (optimization)
        if not self.force and not self.db_interface.confirm_exists(info):
            msg = "No Job with status claimed found in DB for this jobdir. Use force option to bypass this check"
            raise UserWarning(msg)

    def update_db_for_job(self, job_dir_path):
        storage = self.storage_class(job_dir_path, self.name, **self.storage_kwargs)
        try:
            info = storage.load_json(DEFAULT_JOB_INFO_FILENAME)
        except FileNotFoundError as err:
            msg = "No file {} found for job dir {}.".format(DEFAULT_JOB_INFO_FILENAME, storage)
            raise UserWarning(msg)
        except json.decoder.JSONDecodeError as err:
            msg = "JSON Error for job info file in job dir {}.".format(DEFAULT_JOB_INFO_FILENAME, storage)
            raise UserWarning(msg)
        try:
            self._verify_job(info, storage)
            results = self._parse_results_for_completed_job(info, storage)
            self.db_interface.release_job(info, errors=False, force=self.force)
            new_calc_keys = self._add_results_to_db(info, results)
            logger.info("Calcs added to DB: {0}".format(new_calc_keys))
            storage.archive()
            return storage.job_path
        except UserWarning:
            raise
        except FileNotFoundError as err:
            filename = err.filename.lstrip(storage.job_path)
            raise UserWarning("No file {} found for job dir {}.".format(os.path.basename(filename), storage.job_path))
        except dbi.JobNotFoundError as err:
            err.msg += " for dir {}".format(os.path.join(job_dir_path, DEFAULT_JOB_INFO_FILENAME))
            logging.error(err.msg)
            raise UserWarning(err.msg)
        except Exception:
            logger.exception("Failed loading from job dir {0} with loader {1}".format(storage, self.loader_module))
            traceback_string = traceback.format_exc(None)
            storage.dump_file(EXCEPTION_FILENAME, traceback_string)
            storage.mark_error()
            try:
                if info:
                    self.db_interface.release_job(info, errors=True, force=self.force)
            except:
                logger.exception("Unable to release parent after error")
            raise UserWarning('Parse Error, check log file: jobdirparser.log')

    def _parse_results_for_completed_job(self, info, storage):
        with tempfile.TemporaryDirectory() as tempdir:
            storage.copy_to_path(tempdir)
            results = self.loader_module.load_calc_list(tempdir, info)
            if results is None:
                raise Exception("load_calc_list return None")
        return results

    def _add_results_to_db(self, info, results):
        append_num = len(results) > 1
        new_calc_keys = []
        for i, result in enumerate(results):
            try:
                if append_num:
                    suffix = str(i)
                else:
                    suffix = None
                newkey = self.db_interface.add_result(result, info, suffix, self.force)
                new_calc_keys.append(newkey)
            except Exception as e:
                logger.error("Error: {0} when adding result from {1}.".format(e, info))
                raise

        return new_calc_keys

