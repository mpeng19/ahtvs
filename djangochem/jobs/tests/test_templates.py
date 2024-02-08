import unittest
import os

from jinja2 import Template
from jinja2 import Environment, FileSystemLoader, StrictUndefined, DictLoader, UndefinedError

from munch import Munch

MY_DIR = os.path.dirname(__file__)
TEMPLATE_DIR = os.path.join(MY_DIR, "config")
TEMPLATE_FILENAME = "template.txt"

class TestTemplates(unittest.TestCase):
    def setUp(self):
        meta = Munch(smiles="c1ccccc1")
        self.calc = Munch(meta_data=meta)
        atoms = [dict(element='H', x=0, y=5.0, z=6),
                 dict(element='H', x=1, y=6.0, z=6),
                 dict(element='H', x=2, y=7.4, z=6),
                 dict(element='H', x=3, y=8.1, z=6),
                 dict(element='H', x=4, y=9.2, z=6)]

        self.calc.coord_list = atoms

    def test_smiles(self):
        t = Template("{{calc.meta_data.smiles}}")
        result = t.render(calc=self.calc)
        self.assertEquals(result, u"c1ccccc1")

    def test_xyz_output(self):
        expect = """begin atoms
                    H 0.000 5.000 6.000
                    H 1.000 6.000 6.000
                    H 2.000 7.400 6.000
                    H 3.000 8.100 6.000
                    H 4.000 9.200 6.000
                    end atoms"""
        # strip leading spaces from expect
        expect = "\n".join([e.strip() for e in expect.split("\n")])
        env = Environment(loader=FileSystemLoader(TEMPLATE_DIR))
        t = env.get_template(TEMPLATE_FILENAME)
        self.assertEquals(t.render(calc=self.calc), expect)

    def test_strict(self):
        env = Environment(loader=DictLoader({'test': '{{calc.nonexistant_attribute}}'}), undefined=StrictUndefined)
        t = env.get_template("test")
        with self.assertRaises(UndefinedError):
            t.render(calc=self.calc)

    def test_name_filename(self):
        env = Environment(loader=FileSystemLoader(TEMPLATE_DIR))
        t = env.get_template(TEMPLATE_FILENAME)
        self.assertEquals(t.name, TEMPLATE_FILENAME)
        self.assertEquals(t.filename, os.path.join(TEMPLATE_DIR, TEMPLATE_FILENAME))
