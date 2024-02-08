import argparse
import sys

from blocks.block import Block
from django.core.management.base import BaseCommand
from molgen.utils import put_block


class Command(BaseCommand):
    help = ('Add new molecules from SMILES to the database.  '
            'Example: ./manage.py addsmiles test_project molgen/test/molecules.smiles --tag=test')

    def add_arguments(self, parser):
        parser.add_argument('project',
                            help="name of the project")
        parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                            default=sys.stdin,
                            help="file with one SMILES per line or if not specified will be read from keyboard")
        parser.add_argument('-t', '--tag', nargs='+', dest='tags', default=[],
                            help="tags for the new mol objects")

    def handle(self, *args, **options):
        self.main(**options)

    def main(self,
             project,
             infile,
             tags,
             **kwargs):
        self.stdout.write("created, smiles, group, tags")

        for smiles in infile:
            block = Block(smiles=smiles)
            created = put_block(block, project, tags)
            self.stdout.write("{}, {}, {}, {}".format(created, block.smiles(show_markers=False, smiles_stereo=False),
                                                      project, tags))
