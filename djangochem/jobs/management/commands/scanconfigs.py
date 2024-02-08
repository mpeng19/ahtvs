import json
import os
import importlib
import logging
import sys
import traceback
from collections import defaultdict

from django.core.management.base import BaseCommand
import jinja2
from jinja2 import Environment, FileSystemLoader, Undefined
from jinja2.utils import missing, object_type_repr
from jinja2._compat import string_types

from jobs.models import JobConfig
from jobs.dbinterface.dbinterface import WORKER_NAME_KEY, PROJECT_KEY, JOB_KEY


EXPECT_CONFIG_KEYS = ["parent_class_name",
                      "job_template_filename",
                      "loader_module",
                      "name"]

TYPICAL_JOBSPEC = {WORKER_NAME_KEY: 'test_config',
                   JOB_KEY: 1,
                   'priority': 5,
                   'inchikey': 'UHOVQNZJYSORNB-UHFFFAOYSA-N',
                   'smiles': 'c1ccccc1',
                   'details': {'nickname': 'benzene'},
                   PROJECT_KEY: 'test_project',
                   'coords': [{'element': 'C', 'x': 0.000, 'y': 1.396, 'z': 0.000},
                              {'element': 'C', 'x': 1.209, 'y': 0.698, 'z': 0.000},
                              {'element': 'C', 'x': 1.209, 'y': -0.698, 'z': 0.000},
                              {'element': 'C', 'x': 0.000, 'y': -1.396, 'z': 0.000},
                              {'element': 'C', 'x': -1.209, 'y': -0.698, 'z': 0.000},
                              {'element': 'C', 'x': -1.209, 'y': 0.698, 'z':  0.000},
                              {'element': 'H', 'x': 0.000, 'y': 2.479, 'z': 0.000},
                              {'element': 'H', 'x': 2.147, 'y': 1.240, 'z': 0.000},
                              {'element': 'H', 'x': 2.147, 'y': -1.240, 'z': 0.000},
                              {'element': 'H', 'x': 0.000, 'y': -2.479, 'z': 0.000},
                              {'element': 'H', 'x': -2.147, 'y': -1.240, 'z': 0.000},
                              {'element': 'H', 'x': -2.147, 'y': 1.240, 'z':  0.000}],
                   'charge': 0
                   }

CONFIG_ERROR_HEADER = "ERROR -- Config: {}:\n"


class Command(BaseCommand):
    help = 'parse completed jobs into the database'

    def add_arguments(self, parser):
        parser.add_argument('configs_path',
                            type=str,
                            nargs='?',
                            default="../chemconfigs/",
                            help=("directory to scan for json configs. "
                                  "chemconfigs dirs can be nested or not."))
        parser.add_argument('-u', '--update',
                            action='store_true',
                            default=False,
                            help="Update the Database")
        parser.add_argument('--reraise',
                            action='store_true',
                            default=False,
                            help="Reraise exceptions during validation.  Can be useful for debugging templates.")
        parser.add_argument('-d', '--showdetails',
                            action='store_true',
                            default=False,
                            help=("Normally jobspec details are expected to be available in a given Job model.  "
                                  "If you want to highlight these anyway, use this flag to show them as errors"))
        parser.add_argument('-k', '--showkeys',
                            action='store_true',
                            default=False,
                            help=("print all the keys found in all the chemconfig json files "
                                  "(developer tool to help write JobConfig model)"))
    def handle(self, *args, **options):
        self.main(**options)

    def _build_node_path(self, node):
        if hasattr(node, 'name'):
            return node.name
        if hasattr(node, 'node'):
            return self._build_node_path(node.node) + '.' + node.attr
        else:
            return node.attr

    def make_logging_undefined(self):
        """
        taken and modified from jinja2 code
        see: https://github.com/pallets/jinja/blob/master/jinja2/runtime.py
        """
        base = Undefined

        def _log_message(undef):
            if undef._undefined_hint is None:
                if undef._undefined_obj is missing:
                    hint = '{} is undefined'.format(undef._undefined_name)
                elif not isinstance(undef._undefined_name, string_types):
                    hint = '{} has no element {}'.format(
                        object_type_repr(undef._undefined_obj),
                        undef._undefined_name)
                else:
                    hint = '{} has no attribute {}'.format(
                        object_type_repr(undef._undefined_obj), undef._undefined_name)
            else:
                hint = undef._undefined_hint
            attr_nodes = [n for n in self.current_ast.find_all(jinja2.nodes.Getattr)]
            for n in attr_nodes:
                if undef._undefined_name == n.attr:
                    hint += " ({} on line {})".format(self._build_node_path(n), n.lineno)
            name_nodes = [n for n in self.current_ast.find_all(jinja2.nodes.Name)]
            for n in name_nodes:
                if undef._undefined_name == n.name:
                    hint += " (line {})".format(n.lineno)
            self.template_errors.append('Template variable: {}'.format(hint))

        class LoggingUndefined(base):

            def _fail_with_undefined_error(self, *args, **kwargs):
                return base._fail_with_undefined_error(self, *args, **kwargs)

            def __str__(self):
                rv = base.__str__(self)
                _log_message(self)
                return rv

            def __iter__(self):
                rv = base.__iter__(self)
                _log_message(self)
                return rv

            def __bool__(self):
                rv = base.__bool__(self)
                _log_message(self)
                return rv

        return LoggingUndefined

    def _validate_template(self, tname, dirpath, fake_jobspec):
        tpath = os.path.join(dirpath, tname)
        self.template_errors = []
        self.current_template = tname
        try:
            self.current_ast = self.template_env.parse(open(tpath).read())
            self.template_env.get_template(tname).render(jobspec=fake_jobspec)
        except jinja2.exceptions.UndefinedError as err:
            for line in traceback.format_exc(-1).split('\n')[1:]:
                    self.template_errors.append(line)

            if self.raise_errors:
                raise
        except Exception as err:
            self.template_errors.append("Exception during validation (use --reraise to see) : " + str(err))
            if self.raise_errors:
                raise
        if self.template_errors:
            self.errors.append("Template File: " + os.path.join(dirpath, tname))
            for e in self.template_errors:
                self.errors.append("    " + e)

    def _validate_config_contents(self, config, path):
        job_filename = config["job_template_filename"]
        template_filenames = config.get("extra_template_filenames", [])

        dirpath = os.path.dirname(path)

        self.template_env = Environment(loader=FileSystemLoader(dirpath),
                                        undefined=self.make_logging_undefined())
        if job_filename:
            template_name_list = list(template_filenames) + [job_filename]
        else:
            template_name_list = template_filenames

        prep_module_name = config.get("jobspec_prep")

        if prep_module_name is not None:
            module_path = os.path.join(dirpath, prep_module_name + ".py")
            spec = importlib.util.spec_from_file_location(prep_module_name, module_path)
            prep_module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(prep_module)
            prep_func = prep_module.update_jobspec
        else:
            prep_func = None

        fake_jobspec_path = os.path.join(dirpath, 'test', 'fixtures',
                                         'sample_job_spec.json')
        if os.path.exists(fake_jobspec_path):
            fake_jobspec = json.load(open(fake_jobspec_path))
        else:
            fake_jobspec = TYPICAL_JOBSPEC

        if prep_func:
            try:
                prep_func(fake_jobspec)
            except Exception:
                # import pdb; pdb.set_trace()
                # exc_type, exc_obj, exc_tb = sys.exc_info()
                # fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                msg = "Error running prep func (use --reraise to see stack) :"
                self.errors.append(msg)
                for line in traceback.format_exc(-1).split('\n')[1:]:
                    self.errors.append("    " + line)
                if self.raise_errors:
                    raise
        if not self.showdetails:
            fake_jobspec['details'] = defaultdict(str, fake_jobspec['details'])
        for tname in template_name_list:
            self._validate_template(tname, dirpath, fake_jobspec)

        if "default_details_filename" in config:
            if not os.path.isfile(os.path.join(dirpath, config["default_details_filename"])):
                msg = "No default details file found with name {}"
                self.errors.append(msg.format(config["default_details_filename"]))

    def _validate_config(self, config, path):
        self.errors = []
        config_key_set = set(config.keys())
        self.all_config_keys |= config_key_set
        missing_keys = set(EXPECT_CONFIG_KEYS) - config_key_set
        if missing_keys:
            self.errors.append("Missing Keys")
            for k in missing_keys:
                self.errors.append("  - " + k)

        self._validate_config_contents(config, path)

        if self.errors:
            self.stdout.write(CONFIG_ERROR_HEADER.format(path))
            for err in self.errors:
                self.stdout.write("      " + err)
        return self.errors == []

    def _import_configfile(self, path, dryrun=False):
        with open(path) as fp:
            try:
                config = json.load(fp)
                if self._validate_config(config, path):
                    name = config['name']
                    pcn = config['parent_class_name']
                    if not dryrun:
                        jc, created = JobConfig.objects.get_or_create(name=name)
                        if created:
                            self.stdout.write("Adding: {} from path {}".format(name, path))
                            self.created_count += 1
                        else:
                            self.stdout.write("Exists: {} from path {}".format(name, path))
                        jc.parent_class_name = pcn
                        jc.configpath = path
                        jc.save()

                    else:
                        self.stdout.write("Found (dryrun): {} at path {}".format(name, path))
                    self.valid_count += 1
            except json.decoder.JSONDecodeError as err:
                self.stdout.write(CONFIG_ERROR_HEADER.format(path))
                self.stdout.write("  {}".format(err))


    def main(self,
             configs_path,
             update,
             reraise,
             showdetails,
             showkeys,
             **kwargs):
        self.raise_errors = reraise
        self.showdetails = showdetails
        self.created_count = 0
        self.valid_count = 0
        self.all_config_keys = set()
        try_count = 0
        if os.path.isfile(configs_path):
            if configs_path.endswith('config.json'):
                try_count += 1
                self._import_configfile(configs_path, not update)
            else:
                self.stdout.write("filename {} does not end with config.json".format(os.path.basename(configs_path)))
        else:
            for dirpath, dirnames, filenames in os.walk(configs_path):
                for filename in filenames:
                    if filename.endswith('config.json'):
                        path = os.path.join(dirpath, filename)
                        try_count += 1
                        self._import_configfile(path, not update)
        self.stdout.write("Config Count: {}. Valid Count: {}".format(try_count, self.valid_count))
        if update:
            self.stdout.write("Added {} new JobConfigs".format(self.created_count))

        if showkeys:
            self.stdout.write("All keys from chemconfig jsons: ")
            for key in sorted(list(self.all_config_keys)):
                self.stdout.write("    " + key)
