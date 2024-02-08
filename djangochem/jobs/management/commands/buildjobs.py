import json
import os

from tabulate import tabulate

from jobs.jobdirbuilder import JobDirBuilder
from jobs.dbinterface import postgresinterface
from jobs.models import JobConfig, Job

from django.core.management.base import BaseCommand


class Command(BaseCommand):
    help = 'create job directories from Jobs in database'

    def add_arguments(self, parser):
        parser.add_argument('project_name', type=str, default='default')
        parser.add_argument('inbox_path', type=str, default=None)
        parser.add_argument('-c', '--config', dest='config_name', type=str, default=None)
        parser.add_argument('-l', '--limit', type=int, default=None)
        parser.add_argument('-b', '--batchsize', type=int, default=1,
                            help=("build batch jobs with N members each."
                                  "Note that for batch template files, the jobspec is a dictionary with"
                                  "only key-values that are shared amonst all batch members."))
        parser.add_argument('-q', '--query_only', action='store_true', default=False)

    def handle(self, *args, **options):
        self.main(**options)

    def format_jobspec_list(self, jobspec_list):
        """
        take a job list and create a formatted table for printing
        """
        headers = jobspec_list[0].keys()
        headers = sorted(headers)
        headers.pop(headers.index('job_key'))
        headers.insert(0, 'job_key')
        jobspec_list = sorted(jobspec_list, key=lambda k: k['job_key'])
        table = [[row[key] for key in headers] for row in jobspec_list]
        return tabulate(table, headers=headers)

    def main(self, inbox_path, project_name, config_name=None, limit=None, batchsize=1, query_only=False, **kwargs):

        if config_name:
            try:
                jc = JobConfig.objects.get(name=config_name)
            except JobConfig.DoesNotExist:
                self.stdout.write("No JobConfig with name {} was found".format(config_name))
                return
            config_path = jc.configpath
        else:
            query = Job.objects.filter(status__in=['', 'error'], group__name=project_name)
            all_job_configs = query.values_list('config', flat=True).distinct()
            count = len(all_job_configs)
            job_configs = JobConfig.objects.filter(pk__in=all_job_configs)
            if count > 1:
                msg = "Jobs are available with multiple types: {}.  Please specify one with -c"
                self.stdout.write(msg.format([c.name for c in job_configs]))
                return
            elif count < 1:
                self.stdout.write("No jobs found")
                return
            else:
                jc = job_configs[0]
                self.stdout.write("No job config specified.  Proceeding with only option found: {}".format(jc.name))
                config_path = jc.configpath

        config = json.load(open(config_path))
        config_path = os.path.dirname(config_path)
        config_name = config['name']
        db_interface = postgresinterface.PostgresInterface()
        dirbuilder = JobDirBuilder(name=config_name,
                                   project=project_name,
                                   config_path=config_path,
                                   db_interface=db_interface,
                                   job_filename=config["job_template_filename"],
                                   template_filenames=config.get("extra_template_filenames", []),
                                   storage_kwargs={"inbox_job_dir": inbox_path},
                                   jobspec_prep_module_name=config.get("jobspec_prep"),
                                   batch_filename=config.get("batch_template_filename")
                                   )
        if query_only:
            jobspec_list = dirbuilder.get_jobspecs()
            if jobspec_list == []:
                self.stdout.write("No jobs found")
            else:
                table_str = self.format_jobspec_list(jobspec_list)
                self.stdout.write(table_str)
        else:
            count = 0
            job_dirs = dirbuilder.build_job_dirs(limit=limit,
                                                 batch_size=batchsize,
                                                 build_spec=None)
            for job_dir in job_dirs:
                self.stdout.write("created: {}".format(job_dir['storage'].job_path))
                count += 1
            self.stdout.write("Total job directories built: {}".format(count))
