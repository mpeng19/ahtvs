import json
import os

from django.core.management.base import BaseCommand

from jobs.jobdirbuilder import JobDirBuilder
from jobs.models import JobConfig


class Command(BaseCommand):
    help = 'create job directory from job spec'

    def add_arguments(self, parser):
        parser.add_argument('-o', '--output-dir', type=str, default=None,
                            required=True)
        group = parser.add_mutually_exclusive_group(required=True)
        group.add_argument('-c', '--config', dest='config_name', type=str)
        group.add_argument('-f', '--config-file', dest='config_file', type=str)
        parser.add_argument('-s', '--job-spec', dest='jobspec_path', type=str,
                            required=True)

    def handle(self, *args, **options):
        self.main(**options)

    def main(self, output_dir=None, config_name=None, config_file=None,
             jobspec_path=None, **kwargs):
        jobspec = json.load(open(jobspec_path))
        if config_name:
            try:
                jc = JobConfig.objects.get(name=config_name)
                config_path = jc.configpath
            except JobConfig.DoesNotExist:
                self.stdout.write("No JobConfig with name {} was found".format(
                    config_name))
                exit(1)
        elif config_file:
            config_path = config_file
        config = json.load(open(config_path))
        config_path = os.path.dirname(config_path)
        config_name = config['name']
        builder = JobDirBuilder(
            name=config_name,
            config_path=config_path,
            job_filename=config["job_template_filename"],
            template_filenames=config.get("extra_template_filenames", []),
            storage_kwargs={"inbox_job_dir": output_dir},
            jobspec_prep_module_name=config.get("jobspec_prep"),
            batch_filename=config.get("batch_template_filename")
        )
        [job_dir] = builder.build_job_dirs_from_jobspecs(jobspecs=[jobspec])
        self.stdout.write("created: {}".format(job_dir['storage'].job_path))
