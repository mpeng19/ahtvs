import logging
import os
import datetime
import importlib
import time

from jinja2 import Environment, FileSystemLoader, StrictUndefined

from .storage.filestorage import WriteDir
from .dbinterface import dbinterface as dbi

DEFAULT_PRIORITY = 5
DEFAULT_JOB_INFO_FILENAME = "job_info.json"


DEFAULT_TEMPLATE_PATH = "templates"
DEFAULT_PERIOD = 10
DEFAULT_LOG = "worker_error.log"
INBOX_DIR_COUNT_LIMIT = 10
INBOX_DIR_COUNT_TRIGGER = 0

DEFAULT_MIN_PRIORITY = 0
DEFAULT_MAX_PRIORITY = 10

# PRIORITY_RANGE is the max possible range of priorities, tied to the priority system used in the database.
PRIORITY_RANGE = (0, 9)
TEMP_FILE_PREFIX = "temp_"

EXCEPTION_FILENAME = "job_loader_exception.txt"


def common_dict(dict_list):
    """ return dictionary that contains all key value pairs common to *all* the dictionaries in the list"""
    first = dict_list[0]
    common = {}
    for k, v in first.items():
        if type(v) is dict:
            subdicts = [d[k] for d in dict_list if k in d]
            if len(subdicts) == len(dict_list):
                common_sub = common_dict(subdicts)
                common[k] = common_sub
        else:
            matches = [d[k] == v for d in dict_list if k in d]
            if len(matches) == len(dict_list) and all(matches):
                common[k] = v
    return common

class ConfigModuleMissingException(Exception): pass

def import_module_from_config(config_path=None, module_name=None):
    module_path = os.path.join(config_path, module_name + ".py")
    if not os.path.exists(module_path):
        raise ConfigModuleMissingException()
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

class JobDirBuilder(object):
    """ Creates job dirs from database objects requesting work """
    def __init__(self,
                 name,
                 project=None,
                 config_path=None,
                 db_interface=None,
                 min_priority=DEFAULT_MIN_PRIORITY,
                 max_priority=DEFAULT_MAX_PRIORITY,
                 template_filenames=[],
                 job_filename=None,
                 jobspec_prep_module_name=None,
                 batch_filename=None,
                 storage_class=WriteDir,
                 storage_kwargs={}):
        """
        args:
        jobspec_prep_module_name: the name of a module with a function prep_jobspec that accepts
        a parent model and returns a dict that provides extra key:value pairs to be added to the
        jobspec details key
        """
        self.name = name
        self.project = project
        self.config_path = config_path
        self.db_interface = db_interface
        self.min_priority = min_priority
        self.max_priority = max_priority
        self.template_filenames = template_filenames
        self.job_filename = job_filename
        self.batch_filename = batch_filename
        self.storage_class = storage_class
        self.storage_kwargs = storage_kwargs

        try:
            self.job_dir_builder_module = import_module_from_config(
                config_path, module_name='job_dir_builder')
        except ConfigModuleMissingException:
            self.job_dir_builder_module = None

        #### MOVE THIS DOWN!
        self.template_env = Environment(loader=FileSystemLoader(config_path),
                                        undefined=StrictUndefined)
        if batch_filename:
            if os.path.isfile(os.path.join(config_path, batch_filename)):
                self.batch_enabled = True
                self.batch_template = self.template_env.get_template(
                    batch_filename)
                if job_filename is None:
                    raise Exception("The job_filename must be specified"
                                    " to allow for batch mode")
            else:
                msg = "batch filename {} specified but not found in dir {}"
                raise Exception(msg.format(batch_filename, config_path))
        else:
            self.batch_enabled = False


        if job_filename:
            all_template_files = template_filenames + [job_filename]
        else:
            all_template_files = template_filenames
        self.templates = [self.template_env.get_template(tname)
                          for tname in all_template_files]

        if jobspec_prep_module_name is not None:
            self.jobspec_prepper = import_module_from_config(
                config_path, jobspec_prep_module_name)
        else:
            self.jobspec_prepper = None

    def build_job_dirs(self, build_spec=None, limit=None, batch_size=1):
        jobspecs = self.get_jobspecs(limit=limit)
        claimed_jobspecs = self._claim_jobspecs(jobspecs=jobspecs)
        return self.build_job_dirs_from_jobspecs(jobspecs=claimed_jobspecs,
                                              build_spec=build_spec,
                                              batch_size=batch_size)

    def build_job_dirs_from_jobspecs(self, jobspecs=None, build_spec=None,
                                     batch_size=1):
        module_build_job_dirs = getattr(self.job_dir_builder_module,
                                        'build_job_dirs', None)
        if module_build_job_dirs:
            job_dirs = module_build_job_dirs(
                jobspecs=jobspecs,
                build_spec=build_spec,
                gen_storage=self._gen_storage)
        else:
            job_dirs = self._legacy_build_job_dirs(
                jobspecs=jobspecs,
                batch_size=batch_size)
        return job_dirs

    def _gen_storage(self, key=None, extra_storage_kwargs={}):
        merged_storage_kwargs = {**self.storage_kwargs, **extra_storage_kwargs}
        return self.storage_class(key, self.name, **merged_storage_kwargs)

    def _legacy_build_job_dirs(self, jobspecs=[], batch_size=1):
        """
        limit is the max number of jobs to create
        batch_size will batch together n db jobs into a single job dir
        if limit is not a multiple of batch, a final batch with fewer jobs is created
        """
        try:
            job_dirs = []
            if batch_size > 1 and not self.batch_enabled:
                raise Exception("This config is not enabled for batch processing")
            while jobspecs:
                batch = []
                while jobspecs and len(batch) < batch_size:
                    batch.append(jobspecs.pop(0))
                if batch:
                    job_dirs.append(self._build_job_dir(batch))
            return job_dirs
        except:
            for jobspec in jobspecs:
                self.db_interface.clear_job(jobspec)
            logging.error("Clamed jobs cleared. Dirs build:\n".job_dirs)
            raise


    def get_jobspecs(self, limit=None):
        return self.db_interface.get_new_jobs(
            self.name,
            self.project,
            min_priority=self.min_priority,
            max_priority=self.max_priority,
            limit=limit,
            job_to_jobspec=getattr(self.job_dir_builder_module,
                                   'job_to_jobspec', None))
    def _claim_jobspecs(self, jobspecs=[]):
        claimed_jobspecs = []
        for jobspec in jobspecs:
            try:
                self.db_interface.claim_job(jobspec)
                claimed_jobspecs.append(jobspec)
            except dbi.LockError:
                self._log_unclaimable_jobspec(jobspec)
        return claimed_jobspecs

    def _log_unclaimable_jobspec(self, jobspec):
        logging.warning("Skipping jobspec that can't be claimed,"
                       " jobspec was '{}'.".format(jobspec))

    def create_job_from_db(self):
        """create a single job dir from the first job found in db"""
        jobspec = self.get_jobspecs(1)[0]
        try:
            self.db_interface.claim_job(jobspec)
            return self._build_job_dir([jobspec])
        except dbi.LockError:
            self._log_unclaimable_jobspec(jobspec)

    def _make_batch_key_name(self):
        # return "batch_{}_{}".format(datetime.datetime.now().strftime("%Y-%m-%d-%Hh%Mm%S.%f")[:-3],''.join(random.SystemRandom().choice(string.ascii_uppercase + string.digits) for _ in range(4)))
    #    return "batch_{}".format(datetime.datetime.now().strftime("%Y-%m-%d-%Hh%Mm%S.%f")[:-3])
        return "batch_{}".format(str(int(time.time())))

    def _job_key(self, job):
        key = job.get('inchikey')
        if key is None:
            key = job['job_key']
        #return '{}_{}'.format(key,''.join(random.SystemRandom().choice(string.ascii_uppercase + string.digits) for _ in range(4)))
        return '{}'.format(key)

    def _batch_info_dict(self):
        return {dbi.WORKER_NAME_KEY: self.name,
                dbi.PROJECT_KEY: self.project,
                dbi.BATCH_KEY: True}

    def _build_job_dir(self, batch):
        """
        build the job directory and script from a parent document

        if a batch has multiple priorities, the minimum of included jobs used
        """
        if len(batch) == 1:
            key = self._job_key(batch[0])
        else:  # BATCH
            key = self._make_batch_key_name()

        priorities = [jobspec.get('priority', DEFAULT_PRIORITY) for jobspec in batch]
        min_priority = min(priorities)
        storage = self._gen_storage(
            key=key,
            extra_storage_kwargs={'priority': min_priority})
        if len(batch) == 1:
            self._fill_job_dir(storage, batch[0])
        else:  # BATCH
            sub_dir_list = []
            for i, jobspec in enumerate(batch):
                key = self._job_key(jobspec)
                subdirname = "{}_{}".format(i, key)
                sub_dir_list.append(subdirname)
                self._fill_job_dir(storage, jobspec, subdirname)

            storage.dump_json(DEFAULT_JOB_INFO_FILENAME, self._batch_info_dict())
            batchfile_string = self.batch_template.render(jobspec=common_dict(batch),
                                                          sub_dir_list=sub_dir_list,
                                                          job_filename=self.job_filename)
            # write the batch file shell script to the job filename in the root dir
            storage.dump_file(self.job_filename, batchfile_string)
        logging.debug("New job dir created: {0}".format(storage))
        built_job_dir = {'storage': storage}
        return built_job_dir

    def _prepend_if_needed(self, dirname, filename):
        if dirname is None:
            return filename
        else:
            return os.path.join(dirname, filename)

    def _fill_job_dir(self, storage, jobspec, subdirname=None):
        """
        write jobspec info file
        """
        # create jobspec info file

        if subdirname:
            storage.make_subdir(subdirname)

        info_json_path = self._prepend_if_needed(subdirname, DEFAULT_JOB_INFO_FILENAME)
        if self.jobspec_prepper is not None:
            self.jobspec_prepper.update_jobspec(jobspec)

        storage.dump_json(info_json_path, jobspec)

        # build template files
        for template in self.templates:
            rendered_string = template.render(jobspec=jobspec)
            file_path = self._prepend_if_needed(subdirname, template.name)
            storage.dump_file(file_path, rendered_string)
