import json
import os

from tabulate import tabulate

from jobs.jobdirparser import JobDirParser
from jobs.dbinterface import postgresinterface
from jobs.models import JobConfig, Job

from django.core.management.base import BaseCommand

ARCHIVE_DIRNAME = 'archive'
ERROR_DIRNAME = 'error'


class Command(BaseCommand):
    help = 'parse completed jobs into the database'

    def add_arguments(self, parser):
        parser.add_argument('project_name', type=str, default='default')
        parser.add_argument('completed_path', type=str, default=None)
        parser.add_argument('-c', '--config', dest='config_name', type=str, default=None)
        parser.add_argument('-r', '--root_path', type=str, default=None,
                            help="Optional paths to archive and error job dirs. If given, two dirs will be created/used under this path"
                            "'archive' and 'error'.  The default behavior is to leave jobdirs where they were found")
        parser.add_argument('-l', '--limit', type=int, default=None)
        parser.add_argument('-q', '--query_only', action='store_true', default=False)
        parser.add_argument('-f', '--force', action='store_true', default=False,
                            help=("Force parse of job regardless of Job.status. "
                                  "Note, that if you don't pass a config name with -c, force may not work because"
                                  "jobconfig names are discovered by search for open Jobs."))
        parser.add_argument('--lockfile', default=None,
                            help=("""Ignore dirs which contain a file with this
                                  name. The intention of this parameter is to 
                                  avoid errors that come from parsing dirs which
                                  have not been fully copied."""))

    def handle(self, *args, **options):
        self.main(**options)

    def main(self,
             project_name,
             completed_path,
             config_name,
             root_path,
             limit=None,
             query_only=False,
             force=False,
             lockfile=None,
             **kwargs):

        if config_name:
            try:
                job_configs = JobConfig.objects.filter(name=config_name)
            except JobConfig.DoesNotExist:
                self.stdout.write("No JobConfig with name {} was found".format(config_name))
                return
        else:
            query = Job.objects.filter(status__in=['claimed'], group__name=project_name)
            all_job_configs = query.values_list('config', flat=True).distinct()
            count = len(all_job_configs)
            job_configs = JobConfig.objects.filter(pk__in=all_job_configs)
            if count < 1:
                self.stdout.write("No jobs found")
                return

        total_count = 0
        for jc in job_configs:
            self.stdout.write("For JobConfig: {}".format(jc.name))
            config_path = jc.configpath

            config = json.load(open(config_path))
            config_dir_path = os.path.dirname(config_path)
            config_name = config['name']
            loader_module = config['loader_module']

            if root_path is not None:
                root_path = os.path.expanduser(root_path)
                archive_path = os.path.join(root_path, ARCHIVE_DIRNAME)
                error_path = os.path.join(root_path, ERROR_DIRNAME)
            else:
                archive_path = None
                error_path = None
            db_interface = postgresinterface.PostgresInterface()
            dirparser = JobDirParser(name=config_name,
                                     project=project_name,
                                     loader_path=config_dir_path,
                                     loader_module_name=loader_module,
                                     db_interface=db_interface,
                                     force=force,
                                     storage_kwargs={'error_path': error_path,
                                                     'archive_path': archive_path})
            if query_only:
                job_list = dirparser.get_completed_job_list()
                if job_list == []:
                    self.stdout.write("No jobs found")
                else:
                    keys = job_list[0].keys()
                    keys = sorted(keys)
                    job_list = sorted(job_list)
                    self.stdout.write(tabulate(job_list, headers='keys'))
            else:
                count = 0
                # first try to treat this as a jobdir iteself
                try:
                    dirparser.update_db_for_job(completed_path)
                    count += 1
                except UserWarning:
                    pass

                # if no jobinfo file is found, assume it's a completed dir
                if count == 0:
                    for result, error in dirparser.iter_completed_jobs(
                        completed_path, limit, lockfile=lockfile):
                        if error:
                            self.stderr.write(str(error))
                        else:
                            count += 1
                            self.stdout.write("parsed " + result)

                self.stdout.write("{} Jobdirs successfully parsed for {}".format(count, jc.name))
                total_count += count
        if not query_only:
            self.stdout.write("---\n{} Jobdirs successfully parsed across all configs".format(total_count))

