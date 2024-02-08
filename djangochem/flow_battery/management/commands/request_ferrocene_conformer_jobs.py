import argparse
import json
import re

from pgmols.models import Geom
from jobs.management.commands.requestjobs import Command as RequestJobsCommand

FILE_URI_RE = 'file:(.+)'
DB_GEOM_URI_RE = 'db:geom:id:(.+)'

VALID_URIS_HELP = ("Fragment sources must be specified as "
                   "(a) a path to an xyz file: (file:/path/to/my_fragment.xyz)"
                   " or (b) a database geom id: (db:geom:id:1234).")

def fragment_source_arg_handler(string):
    if re.match(FILE_URI_RE, string) or re.match(DB_GEOM_URI_RE, string):
        return string
    raise argparse.ArgumentTypeError("Invalid fragment source '{input}'."
                                     "{valid_uris_help}".format(
                                         input=string,
                                         valid_uris_help=VALID_URIS_HELP))

class Command(RequestJobsCommand):
    help = 'Request ferrocene conformer jobs.'

    def add_arguments(self, parser):
        # should take args for 2 frag sources.
        # frag sources can specified as some kind of uri, e.g 'file:asdfasf', or
        # 'db:geom:id:someId'.
        for i in [1, 2]:
            parser.add_argument( '--frag' + str(i),
                                help=("Source for fragment geometry {frag_num}."
                                      " {valid_uris_help}").format(
                                          frag_num=1,
                                          valid_uris_help=VALID_URIS_HELP),
                                required=True)
        # Add args for normal requestjobs.
        RequestJobsCommand.add_arguments(self, parser)

    def handle(self, *args, **options):
        self.main(**options)

    def main(self, frag1=None, frag2=None, **kwargs):
        fragment_xyzs = [self.resolve_fragment_uri_as_xyz(frag_uri)
                         for frag_uri in [frag1, frag2]]
        self.update_details(kwargs, fragment_xyzs)
        RequestJobsCommand.main(self, **kwargs)

    def update_details(self, kwargs, fragment_xyzs):
        details = json.loads(kwargs.get('details', '{}'))
        fragment_details = {
            'fragments': [
                {
                    'key': 'frag1',
                    'xyz': fragment_xyzs[0],
                },
                {
                    'key': 'frag2',
                    'xyz': fragment_xyzs[1],
                }
            ]
        }
        details.update(fragment_details)
        kwargs['details'] = json.dumps(details)

    def resolve_fragment_uri_as_xyz(self, frag_uri):
        try:
            file_match = re.match(FILE_URI_RE, frag_uri)
            if file_match:
                return open(file_match.group(1)).read()
            db_geom_match = re.match(DB_GEOM_URI_RE, frag_uri)
            if db_geom_match:
                geom = Geom.objects.get(id=db_geom_match.group(1))
                return geom.as_xyz()
        except Exception as ex:
            raise Exception("Could not load fragment from uri '{frag_uri}':"
                            " {ex}".format(frag_uri=frag_uri, ex=ex))
