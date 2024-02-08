import os

from django.contrib.auth.models import Group
from django.test import TestCase

from jobs.jobdirbuilder import common_dict


class TestJobspecs(TestCase):
    def setUp(self):
        pass

    def test_common(self):
        a = {'a': 1,
             'b': 2,
             'c': 3}

        b = {'a': 1}
        c = {'k': 5,
             'c': 3,
             'b': 2}

        d = {'a': {'b': 1, 'c': 2}}
        e = {'a': {'b': 1, 'c': 3}}

        expect = {'a': 1}
        common = common_dict([a,b])
        self.assertEqual(expect, common)

        expect = {}
        common = common_dict([a,b,c])
        self.assertEqual(expect, common)

        expect = {'b': 2, 'c': 3}
        common = common_dict([c,a])
        self.assertEqual(expect, common)

        expect = {'a': {'b': 1}}
        common = common_dict([d,e])
        self.assertEqual(expect, common)
