import json

from jobs import jobrequester

from tabulate import tabulate

from django.core.management.base import BaseCommand
from jobs.models import JobConfig


class Command(BaseCommand):
    help = 'Request new jobs: create Jobs entries in the database'

    def add_arguments(self, parser):
        parser.add_argument('project_name', type=str,
                            help='the project to work within (root Mols must be members)')
        parser.add_argument('request_config', type=str,
                            help="the chemconfig to request on the parent objects")
        parser.add_argument('-p', '--parent_config', type=str, default=None,
                            help="Limit to models who were created from jobs of this config")
        parser.add_argument('-l', '--limit', type=int, default=None,
                            help="Only create this many new jobs.")
        parser.add_argument('--workbatch', type=str, default=None,
                            help=("Limit job requests to members of this workbatch. "
                                  "New jobs will become members of this workbatch as well."))
        parser.add_argument('-t', '--tag', nargs='+', dest='tags', default=[],
                            help="only run on molecules (or their children) who have these tags")
        parser.add_argument('--clear_errors', action='store_true', default=False,
                            help=("Find jobs of requested config name with error status,"
                                  " and clear the error status to allow re-building."))
        parser.add_argument('--details', type=str, default='{}',
                            help='a json string to become details of the requested jobs')
        parser.add_argument('--allow_dupes', action='store_true', default=False,
                            help="allow the same config to be called again on object (very unusual)")
        parser.add_argument('--force', action='store_true', default=False,
                            help="Force through large orders (over 10000 jobs).  eg. bypass warnings")
        parser.add_argument('-s', '--request-specs', dest='request_spec_paths',
                            nargs='*', help=("List of request spec json files."
                                             " Multiple specs can be provided,"
                                             " with properties in later specs"
                                             " overriding earlier ones."))
        parser.add_argument('--dry_run', action='store_true', default=False,
                            help=("Don't claim jobs, just show how many would"
                                  " be created."))

    def handle(self, *args, **options):
        self.main(**options)

    def format_job_list(self, job_list):
        """
        take a job list and create a formatted table for printing
        """
        headers = ['config', 'status']
        table = [(job.config.name, job.status) for job in job_list]
        return tabulate(table, headers=headers)

    def main(self,
             project_name=None,
             request_config=None,
             parent_config=None,
             limit=None,
             workbatch=None,
             tags=None,
             details=None,
             allow_dupes=None,
             verbosity=None,
             force=None,
             clear_errors=None,
             request_spec_paths=None,
             dry_run=None,
             **kwargs):

        if clear_errors:
            result = jobrequester.clear_error_status(project_name,
                                                     request_config,
                                                     limit=limit)
            self.stdout.write("{} Jobs with status error were cleared".format(result))
        else:
            try:
                jobs = jobrequester.request_jobs(
                    project_name,
                    request_config,
                    parent_config_name=parent_config,
                    limit=limit,
                    workbatch=workbatch,
                    details=json.loads(details),
                    avoid_dupes=not allow_dupes,
                    tags=tags,
                    force=force,
                    request_spec=self.combine_request_specs(
                        request_spec_paths),
                    dry_run=dry_run,
                )
                if dry_run:
                    self.stdout.write("[DRY_RUN] Would have requested "
                                      "{} new jobs".format(len(jobs)))
                else:
                    self.stdout.write("Requested {} new jobs".format(len(jobs)))
                if verbosity > 2:
                    self.stdout.write("Requested Jobs:")
                    self.stdout.write(self.format_job_list(jobs))

            except JobConfig.DoesNotExist:
                msg = ("Requested config '{}' Does not exist.  "
                       "You may need to run scanconfig to add this config to the database.")
                self.stdout.write(msg.format(request_config))

    def combine_request_specs(self, request_spec_paths):
        if request_spec_paths is None:
            return None
        request_specs = [json.load(open(request_spec_path))
                         for request_spec_path in request_spec_paths]
        combined_spec = {}
        for spec in request_specs:
            combined_spec.update(spec)
        return combined_spec

