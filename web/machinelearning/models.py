from django.db import models
from django.utils import timezone
from django.db.models import When, Value, Case
from django.forms.models import model_to_dict

from rdkit.Chem import AllChem
from autograd import numpy as np

from molvote.models import Method, Property, Candidate
from ml import nnautograd


FINGERPRINT_SIZE = 512
DEFAULT_HIDDEN_SIZE = 100
DEFAULT_HIDDEN_COUNT = 1
DEFAULT_SCALE = 0.1
DEFAULT_MOMENTUM = 0.98371
DEFAULT_LEARNRATE = 0.001604
DEFAULT_EPOCHS = 50
DEFAULT_BATCHSIZE = 256
DEFAULT_L2 = 0.0

def smile_to_fp(s, fplength, radius):
    """Computes fingerprints from a SMILES string."""
    m = AllChem.MolFromSmiles(s)
    bitstring = AllChem.GetMorganFingerprintAsBitVect(m, radius, nBits=fplength).ToBitString()
    return np.array(list(bitstring), dtype=int)


def get_props(candidates, prop_name, method_contains="s0"):
    method = Method.objects.get(name__contains=method_contains)
    train_props = Property.objects.filter(candidate__in=candidates,
                                          method=method,
                                          name=prop_name)
    train_props = train_props.select_related("candidate")
    return train_props


def props_to_fingerprints_and_vals(props):
    smiles = [p.candidate.smiles for p in props]
    vals = np.array([p.value for p in props])
    fingerprints = np.array([smile_to_fp(s, fplength=512, radius=4) for s in smiles])
    return fingerprints, vals


def unpack_fingerprint(bits):
    return np.unpackbits(np.fromstring(bits, dtype=np.uint8))


def query_fp_and_vals(qs, prop_name, qc_method, fp_method, remove_low_strength=True):
    qs = qs.annotate(prop_val=Case(When(properties__name=prop_name,
                                        properties__method=qc_method,
                                        then='properties__value'),
                                   default=Value(None))).exclude(prop_val=None)
    qs = qs.annotate(fp=Case(When(fingerprint__method=fp_method,
                                  then='fingerprint__fingerprint'),
                             default=Value(None))).exclude(fp=None)
    if remove_low_strength:
        qs = qs.filter(properties__name='strength',
                       properties__method=qc_method,
                       properties__value__gt=0.00001)
    return qs.values_list('fp', 'prop_val')



class NeuralNet(models.Model):
    weights = models.BinaryField()
    mean = models.FloatField(null=True)
    std = models.FloatField(null=True)
    name = models.CharField(max_length=255, null=True, blank=True)
    created_at = models.DateTimeField(default=timezone.now)

    training_size = models.IntegerField(default=0)
    batch_size = models.IntegerField(default=DEFAULT_BATCHSIZE)
    epochs = models.IntegerField(default=DEFAULT_EPOCHS)
    init_scale = models.FloatField(default=DEFAULT_SCALE)
    momentum = models.FloatField(default=DEFAULT_MOMENTUM)
    learn_rate = models.FloatField(default=DEFAULT_LEARNRATE)
    input_size = models.IntegerField(default=FINGERPRINT_SIZE)
    hidden_size = models.IntegerField(default=DEFAULT_HIDDEN_SIZE)
    hidden_layers = models.IntegerField(default=DEFAULT_HIDDEN_COUNT)
    l2_reg = models.FloatField(default=DEFAULT_L2)
    qc_method = models.ForeignKey(Method, related_name='qc_method', null=True)
    prop_name = models.CharField(max_length=255, null=True, blank=True)
    fp_method = models.ForeignKey(Method, related_name='fp_method', null=True)

    def train(self, candidate_query, qc_method, prop_name, fp_method, progress_func=None, test_data=None):
        self.training_fp_and_vals = query_fp_and_vals(candidate_query,
                                                      prop_name,
                                                      qc_method,
                                                      fp_method)

        if len(self.training_fp_and_vals) == 0:
            raise StandardError("no data found with fingerprints and specified properties")
        # unzip
        fingerprints, vals = zip(*self.training_fp_and_vals)
        self.qc_method = qc_method
        self.fp_method = fp_method
        self.prop_name = prop_name
        self.train_with_raw_fp(fingerprints, vals, progress_func, test_data)

    def train_with_raw_fp(self, raw_fingerprints, values, progress_func=None, test_data=None):
        fingerprint_arrays = np.array([unpack_fingerprint(fp) for fp in raw_fingerprints])
        vals = np.array(values)
        self.num_inputs = fingerprint_arrays.shape[1]

        if test_data:
            test_fps, test_vals = zip(*test_data)
            test_fps = [unpack_fingerprint(fp) for fp in test_fps]
            test_data = zip(test_fps, test_vals)

        param_dict = nnautograd._train(features=fingerprint_arrays,
                                       raw_targets=vals,
                                       h1_size=self.hidden_size,
                                       h1_layer_count=self.hidden_layers,
                                       momentum=self.momentum,
                                       learn_rate=self.learn_rate,
                                       l2_reg=self.l2_reg,
                                       init_scale=self.init_scale,
                                       batch_size=self.batch_size,
                                       num_epochs=self.epochs,
                                       progress_func=progress_func,
                                       test_data=test_data)

        self.raw_weights = param_dict["weights"]
        self.weights = param_dict["weights"].astype(np.float64).tostring()
        self.mean = param_dict["mean"]
        self.std = param_dict["std"]
        self.training_size = vals.size

    def predict(self, fingerprints, assert_fp_type_match=False):
        """return prediction from array of fingerprint objects"""
        if assert_fp_type_match:
            assert(all([fp.name == self.fingerprint_name for fp in fingerprints]))
        raw_fingerprints = [fp.fingerprint for fp in fingerprints]
        return self.predict_with_raw_fp(raw_fingerprints)

    def predict_with_raw_fp(self, raw_fingerprints):
        fingerprint_arrays = np.array([unpack_fingerprint(fp) for fp in raw_fingerprints])
        weights_arr = np.fromstring(self.weights, dtype=np.float64)
        return nnautograd._predict(w_vect=weights_arr,
                                   mean=self.mean,
                                   std=self.std,
                                   features=fingerprint_arrays,
                                   h1_size=self.hidden_size,
                                   h1_layer_count=self.hidden_layers)

    def dump_to_dict(self):
        out_dict = {}
        for f in self._meta.get_fields():
            if hasattr(f, 'attname'):
                out_dict[f.attname] = getattr(self, f.attname)
        out_dict['qc_method'] = self.qc_method.name
        out_dict['fp_method'] = self.fp_method.name
        return out_dict

    def load_from_dict(self, input_dict):
        for key, value in input_dict.items():
            if key not in ['qc_method', 'fp_method']:
                setattr(self, key, value)
        self.qc_method = Method.objects.get(name=input_dict['qc_method'])
        self.fp_method = Method.objects.get(name=input_dict['fp_method'])


class Fingerprint(models.Model):
    fingerprint = models.BinaryField()
    method = models.ForeignKey(Method, null=True)
    candidate = models.ForeignKey(Candidate, null=True)

    def _morgan_method_name(self, radius, fplength):
        return "morgan-{}-r{}".format(fplength, radius)

    def morgan(self, candidate=None, radius=4, fplength=512, method=None):
        if method is None:
            if self.method is None:
                self.method, created = Method.objects.get_or_create(name=self._morgan_method_name(radius, fplength))
        else:
            self.method = method
        if candidate:
            self.candidate = candidate
        self._morgan_smiles(self.candidate.smiles, radius, fplength)

    def _morgan_smiles(self, smiles, radius, fplength):
        if fplength % 8 != 0:
            raise StandardError("fingerprints length must be multiples of 8")
        m = AllChem.MolFromSmiles(smiles)
        bit_string = AllChem.GetMorganFingerprintAsBitVect(m, radius, nBits=fplength).ToBitString()
        arr = np.array(list(bit_string), dtype=np.uint8)
        self._store_array(arr)

    def _store_array(self, array):
        """ compress array of ints or bools to single bits and write to storage"""
        self.fingerprint = np.packbits(array).tostring()

    def array(self):
        """ retrieve stored bits and expand to unambigious array of uint8"""
        return unpack_fingerprint(self.fingerprint)
