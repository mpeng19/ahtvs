import logging

from molgen import blockassembler
from django.core.management.base import BaseCommand
from rdkit import RDLogger


class Command(BaseCommand):
    help = ('Generate new molecules from fragments and add to database.  '
            'Example: ./manage.py molgen test_project molgen/test/fragments_small.json --recipe=test')

    def add_arguments(self, parser):
        parser.add_argument('project',
                            help="name of the project")
        parser.add_argument('inputfile', type=str, default="tests/fragments.json",
                            help='A JSON file containing the initial fragments')
        parser.add_argument(
            '-r', '--recipe', type=str,
            choices=['oled', 'cep', 'laser', 'test', 'ferrocene', 'spiro',
                     'symmetric_terminal_spacer_core','noble_graft'],
            default='oled', help="select the way fragments are combined")
        parser.add_argument('--style', type=str, choices=['fill', 'all', 'sym'], default='sym',
                            help="use react_all or react_fill ontop of react_sym")
        parser.add_argument('-t', '--tag', nargs='+', dest='tags', default=[],
                            help="tags for the new mol objects")
        parser.add_argument("-l",
                            "--log",
                            dest="loglevel",
                            default='WARNING',
                            choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                            help="Set the logging level")
        parser.add_argument("--rdlog",
                            dest="rdloglevel",
                            default='ERROR',
                            choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                            help="Set the logging level for rdkit reactions (def: ERROR)")

    def handle(self, *args, **options):
        self.main(**options)

    def main(self,
             project,
             inputfile,
             recipe,
             style,
             tags,
             loglevel,
             rdloglevel,
             **kwargs):

        RDLogger.logger().setLevel(getattr(RDLogger, rdloglevel))

        logger = logging.getLogger()
        loglevel = getattr(logging, loglevel)
        logger.setLevel(loglevel)
        ch = logging.StreamHandler(self.stdout)
        ch.setLevel(loglevel)

        blockassembler.run(inputfile, project, recipe, self.stdout, tags=tags, style=style)
