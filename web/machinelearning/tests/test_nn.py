import autograd.numpy.random as npr
import numpy.testing as npt
import numpy as np
from django.test import TestCase
from django.db.models import F, Q, When, Value, Case

from rdkit.Chem import AllChem

from molvote.models import Property, CandidateSet, Method, Candidate
from machinelearning.models import NeuralNet, Fingerprint, unpack_fingerprint, query_fp_and_vals




SET_NAME = "experimental_calibration"


class TestNN(TestCase):
    fixtures = ['machinelearning/tests/fixture.json']

    def setUp(self):
        npr.seed(0)
        cset = CandidateSet(name="test")
        cset.save()
        for c in Candidate.objects.filter(sets__name="experimental_calibration")[:20]:
            c.sets.add(cset)
        self.cans = Candidate.objects.filter(sets=cset)

    def test_query(self):
        props = Property.objects.filter(candidate__sets__name=SET_NAME, method__name__contains="s0", name="splitting")
        self.assertEquals(50, props.count())

    def _train_nn(self):
        for c in self.cans:
            f = Fingerprint()
            f.morgan(c, fplength=128)
            f.save()
        qc_method = Method.objects.get(name__contains="s0")
        fp_method = Method.objects.get(name__contains="morgan")
        self.nn = NeuralNet(batch_size=1, epochs=20)
        self.nn.train(candidate_query=self.cans,
                      prop_name="splitting",
                      qc_method=qc_method,
                      fp_method=fp_method)
        self.nn.save()
        return qc_method, fp_method


    def _train_deep_nn(self):
        for c in self.cans:
            f = Fingerprint()
            f.morgan(c, fplength=128)
            f.save()
        qc_method = Method.objects.get(name__contains="s0")
        fp_method = Method.objects.get(name__contains="morgan")
        self.nn = NeuralNet(batch_size=1, hidden_layers=2, epochs=50)
        self.nn.train(candidate_query=self.cans,
                      prop_name="splitting",
                      qc_method=qc_method,
                      fp_method=fp_method)
        self.nn.save()
        return qc_method, fp_method


    def test_predict(self):
        qc_method, fp_method = self._train_nn()
        nn = NeuralNet.objects.all().first()
        fp_and_vals = query_fp_and_vals(self.cans, "splitting", qc_method, fp_method)
        fps, expect = zip(*fp_and_vals)
        npt.assert_array_almost_equal(nn.predict_with_raw_fp(fps), expect, decimal=2)


    def test_deep_predict(self):
        qc_method, fp_method = self._train_deep_nn()
        nn = NeuralNet.objects.all().first()
        fp_and_vals = query_fp_and_vals(self.cans, "splitting", qc_method, fp_method)
        fps, expect = zip(*fp_and_vals)
        npt.assert_array_almost_equal(nn.predict_with_raw_fp(fps), expect, decimal=2)

    # def train_linreg(self):
    #     qc_method = Method.objects.get(name__contains="s0")
    #     fp_method = Method.objects.get(name__contains="morgan")
    #     self.nn = NeuralNet(num_hidden=0)
    #     self.nn.train(candidate_query=self.cans,
    #                   prop_name="splitting",
    #                   qc_method=qc_method,
    #                   fp_method=fp_method)
    #     self.nn.save()
    #     return qc_method, fp_method

    # def test_linreg_predict(self):
    #     qc_method, fp_method = self.train_linreg()
    #     nn = NeuralNet.objects.all().first()
    #     fp_and_vals = query_fp_and_vals(self.cans, "splitting", qc_method, fp_method)
    #     fps, expect = zip(*fp_and_vals)
    #     npt.assert_array_almost_equal(nn.predict_with_raw_fp(fps), expect, decimal=2)

    def test_fingerprints(self):
        """ test the use of annotation in the current setup to get a set of fingerprints and properties
        in a single database call"""

        method = Method.objects.get(name__contains="s0")
        fp_method, created = Method.objects.get_or_create(name="morgan-512-r4")
        fp_method.save()
        for c in list(Candidate.objects.all()):
            fgrprnt, created = Fingerprint.objects.get_or_create(candidate=c, method=fp_method)
            if created:
                fgrprnt.morgan()
                fgrprnt.save()

        qs = Candidate.objects.all()
        qs = qs.annotate(prop_val=Case(When(properties__name="splitting",
                                            properties__method=method,
                                            then='properties__value'),
                                       default=Value(None))).exclude(prop_val=None)
        qs = qs.annotate(fp=Case(When(fingerprint__method__name="morgan-512-r4",
                                      then='fingerprint__fingerprint'),
                                 default=Value(None))).exclude(fp=None)
        data = qs.values('fp', 'prop_val')
        data = dict([(tuple(unpack_fingerprint(c['fp'])), c['prop_val']) for c in data])

        self.assertEquals(len(data), 136)

        for i in range(5):
            test_prop = Property.objects.filter(name="splitting",
                                                method=method).order_by("?").first()
            # if Fingerprint.objects.filter(candidate=test_prop.candidate).count() > 1:
            #     import ipdb; ipdb.set_trace()
            fp = Fingerprint.objects.get(candidate=test_prop.candidate, method=fp_method)
            self.assertEquals(test_prop.value, data[tuple(fp.array())])

    def tearDown(self):
        Fingerprint.objects.all().delete()
        NeuralNet.objects.all().delete()


class TestFingerprint(TestCase):
    def setUp(self):
        npr.seed(0)
        self.can = Candidate()
        self.can.save()
        self.smiles = "[Mg]N1C2=CC=C([Be])C=C2C2=C1C=CC([Be])=C2"
        self.radius = 4
        fplength = 512
        mol = AllChem.MolFromSmiles(self.smiles)
        expect = AllChem.GetMorganFingerprintAsBitVect(mol, self.radius, nBits=fplength).ToBitString()
        self.expect = np.array(list(expect), dtype=np.uint8)

    def test_morgan(self):
        fplength = 512
        fp = Fingerprint(candidate=self.can)
        fp._morgan_smiles(self.smiles, self.radius, fplength)
        self.assertTrue((self.expect == fp.array()).all())
        self.assertEqual(len(fp.fingerprint), fplength/8)
        return fp

    def test_morgan_odd_size(self):
        fplength = 513
        fp = Fingerprint(candidate=self.can)
        with self.assertRaises(StandardError):
            fp._morgan_smiles(self.smiles, self.radius, fplength)

    def test_morgan_retrieve(self):
        fp = self.test_morgan()
        fp.save()
        fp = Fingerprint.objects.all().first()
        self.assertTrue((self.expect == fp.array()).all())

    def tearDown(self):
        Fingerprint.objects.all().delete()


