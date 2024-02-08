import os
import hashlib
import base64
from collections import defaultdict

from django.utils import timezone

from django.db import models
from django.contrib.auth.models import User


class AccessControlled(models.Model):
    class Meta:
        abstract = True
    project = models.CharField(max_length=255, blank=True, default="")
    released = models.BooleanField(default=False)  # sets whether non-super-users can see this set


class Candidate(AccessControlled):
    """ analog to molecule, but we don't want to reuse that class name.  too confusing!"""
    inchi_key = models.CharField(max_length=27)
    smiles = models.CharField(max_length=1000)
    absorption = models.FloatField(null=True)
    splitting = models.FloatField(null=True)
    strength = models.FloatField(null=True)
    rate = models.FloatField(null=True)
    homo = models.FloatField(null=True)
    lumo = models.FloatField(null=True)
    nicknames = models.CharField(max_length=1000, blank=True, default="")  # comma separated list
    calc_time = models.DateTimeField(null=True)
    sascore = models.FloatField(null=True)
    weight = models.FloatField(null=True)
    total_energy = models.FloatField(null=True)
    water_solvation_energy = models.FloatField(null=True)
    sets = models.ManyToManyField('CandidateSet')
    cas = models.CharField(max_length=12, blank=True, null=True)
    donor = models.ForeignKey("Candidate", null=True, related_name="donor_for")
    acceptor = models.ForeignKey("Candidate", null=True, related_name="acceptor_for")
    bridge1 = models.ForeignKey("Candidate", null=True, related_name="bridge1_for")
    bridge2 = models.ForeignKey("Candidate", null=True, related_name="bridge2_for")

    def __unicode__(self):
        return self.inchi_key

    def short_inchi(self):
        return self.inchi_key.split("-")[0]

    def props(self):
        out = defaultdict(dict)
        for p in Property.objects.filter(candidate=self):
            out[p.method.name][p.name] = p.value
        return out


class Method(models.Model):
    THEORY_MM = "mm"
    THEORY_SE = "se"
    THEORY_DFT = "dft"
    THEORY_WF = "wf"
    THEORY_CHOICES = [("mm", "Molecular Mechanics"),
                        ("se", "Semiempirical"),
                        ("dft", "DFT"),
                        ("wf", "Wave Function")]

    LEVEL_GROUNDSTATE = "gs"
    LEVEL_RPA = "rpa"
    LEVEL_TDA = "tda"
    LEVEL_MMFF94 = "mmff94"
    LEVEL_UFF = "uff"
    LEVEL_PM7 = "pm7"

    LEVEL_CHOICES = [("DFT levels", [("gs", "Ground State"),
                                     ("rpa", "RPA"),
                                     ("tda", "TDA")]),
                     ("MM Forcefields", [("mmff94", "MMFF94"),
                                         ("uff", "UFF")]),
                     ("Semiempirical method",[("pm7", "PM7")])
                     ]

    BASIS_631_GS = "631gs"
    BASIS_6311_PGDP = "6311pgdp"
    BASIS_CHOICES = [("631gs", "6-31G*"),
                     ("6311pgdp", "6-311+Gdp")]

    FUNCTIONAL_B3LYP = "b3lyp"
    FUNCTIONAL_M062X = "m062x"
    FUNCTIONAL_WB97XD = "wb97xd"
    FUNCTIONAL_PBE = "pbe"
    FUNCTIONAL_CHOICES = [("b3lyp", "B3LYP"),
                             ("m062x", "M062X"),
                             ("wb97xd", "Omega-B97xD"),
                             ("pbe", "PBE")]

    SOLVATION_MODEL_IEFPCM = "iefpcm"
    SOLVATION_MODEL_CPCM = "cpcm"
    SOLVATION_MODEL_COSMO = "cosmo"
    SOLVATION_MODEL_CHOICES = [("iefpcm", "IEFPCM"),
                               ("cpcm", "CPCM"),
                               ("cosmo", "COSMO")]

    SOLVENT_WATER = "water"
    SOLVENT_TOLUENE = "toluene"
    SOLVENT_CHOICES = [("water", "Water"),
                       ("toluene", "Toluene")]

    SHELL_RESTRICTED_CLOSED = "r"
    SHELL_UNRESTRICTED_OPEN = "u"
    SHELL_RESTRICTED_OPEN = "ro"

    SHELL_CHOICES = [("r", "Restricted"),
                     ("u", "Unrestricted"),
                     ("ro", "Openshell")]

    MULTIPLICITY_CHOICES = [(1, "singlet"),
                            (2, "doublet"),
                            (3, "triplet")]

    geom_theory = models.CharField(max_length=12, blank=True, choices=THEORY_CHOICES)
    geom_level = models.CharField(max_length=12, blank=True, choices=LEVEL_CHOICES)
    geom_basis = models.CharField(max_length=12, blank=True, choices=BASIS_CHOICES)
    geom_functional = models.CharField(max_length=12, blank=True, choices=FUNCTIONAL_CHOICES)
    geom_solvation = models.CharField(max_length=12, blank=True, choices=SOLVATION_MODEL_CHOICES)
    geom_solvent = models.CharField(max_length=12, blank=True, choices=SOLVENT_CHOICES)
    geom_shell = models.CharField(max_length=12, blank=True, choices=SHELL_CHOICES)
    geom_multiplicity = models.IntegerField(null=True, choices=MULTIPLICITY_CHOICES)

    sp_theory = models.CharField(max_length=12, blank=True, choices=THEORY_CHOICES)
    sp_level = models.CharField(max_length=12, blank=True, choices=LEVEL_CHOICES)
    sp_basis = models.CharField(max_length=12, blank=True, choices=BASIS_CHOICES)
    sp_functional = models.CharField(max_length=12, blank=True, choices=FUNCTIONAL_CHOICES)
    sp_solvation = models.CharField(max_length=12, blank=True, choices=SOLVATION_MODEL_CHOICES)
    sp_solvent = models.CharField(max_length=12, blank=True, choices=SOLVENT_CHOICES)
    sp_shell = models.CharField(max_length=12, blank=True, choices=SHELL_CHOICES)
    sp_multiplicity = models.IntegerField(null=True, choices=MULTIPLICITY_CHOICES)

    name = models.CharField(max_length=255, unique=True)

    def __unicode__(self):
        return self.name

class Property(models.Model):
    """A property is a value relating a calculationtype with a candidate structure"""

    candidate = models.ForeignKey(Candidate, related_name="properties")
    method = models.ForeignKey(Method)
    name = models.CharField(max_length=255)
    value = models.FloatField()


class RedoxPair(AccessControlled):
    reduced = models.ForeignKey(Candidate, related_name="reduced_of_redox_pair")
    oxidized = models.ForeignKey(Candidate, related_name="oxidized_of_redox_pair")
    hydrated_oxidized = models.ForeignKey(Candidate, null=True, related_name="hyd_oxidized_of_redox_pair")
    michael_hydration_energy = models.FloatField(null=True)
    log_hyd_constant = models.FloatField(null=True)
    log_hyd_constant_error = models.FloatField(null=True)
    redox_potential = models.FloatField(null=True)
    redox_potential_error = models.FloatField(null=True)
    is_minimum = models.BooleanField(default=False)
    water_solvation_energy = models.FloatField(null=True)
    substituent_count = models.IntegerField(null=True)

    # reduced_water_solvation_energy = models.FloatField(null=True)

    def family_graph(self):
        current_reduced = self.reduced
        while current_reduced.oxidized_of_redox_pair.all():
            pair_to_left = current_reduced.oxidized_of_redox_pair.first()
            current_reduced = pair_to_left.reduced
        edges = []
        nodes = []
        current_reduced.depth = 0
        nodes_to_search = [current_reduced]
        while nodes_to_search:
            red = nodes_to_search.pop()
            try:
                reduced_index = nodes.index(red)
            except ValueError:
                nodes.append(red)
                reduced_index = len(nodes) - 1
            for pair in red.reduced_of_redox_pair.all():
                ox = pair.oxidized
                ox.depth = red.depth + 1
                try:
                    ox_index = nodes.index(ox)
                except ValueError:
                    nodes.append(ox)
                    ox_index = len(nodes) - 1
                edges.append({"source": reduced_index, "target": ox_index})
                nodes_to_search.append(ox)
        # turn list of inchi keys into list of dicts
        nodes = [{"inchi_key":mol.inchi_key, "depth": mol.depth} for mol in nodes]
        return {"nodes":nodes, "edges":edges}

    def _add_child(self, reduced, is_minimum):
        tree = {"inchi_key": reduced.inchi_key, "min": is_minimum}
        children = reduced.reduced_of_redox_pair.all()
        if children:
            tree["children"] = []
            for pair in children:
                ox = pair.oxidized
                tree["children"].append(self._add_child(ox, is_minimum and pair.is_minimum))
        return tree

    def family_tree(self):
        current_reduced = self.reduced
        while current_reduced.oxidized_of_redox_pair.all():
            pair_to_left = current_reduced.oxidized_of_redox_pair.first()
            current_reduced = pair_to_left.reduced
        return self._add_child(current_reduced, is_minimum=True) # root node minumum by default


class RedoxPairOfPair(AccessControlled):
    low_rp_pair = models.ForeignKey(RedoxPair, related_name="low_rp_of_pair")
    high_rp_pair = models.ForeignKey(RedoxPair, related_name="high_rp_of_pair")
    redox_potential = models.FloatField(null=True)
    redox_potential_error = models.FloatField(null=True)


def make_key():
    return base64.urlsafe_b64encode(os.urandom(12))


class CandidateSet(AccessControlled):
    """ analog to tag or molecule set"""
    name = models.CharField(max_length=255)
    #members = models.ManyToManyField(Candidate)
    creation_time = models.DateTimeField(null=True)
    announced = models.BooleanField(default=False)
    creator = models.ForeignKey(User, null=True)
    key = models.CharField(max_length=100, null=True, blank=True, unique=True, default=make_key)

    def __unicode__(self):
        return self.name

    def size(self):
        return Candidate.objects.filter(sets=self).count()

    def announce(self):
        self.announced = True
        self.creation_time = timezone.now()
        self.release()

    def release(self):
        self.released = True
        self.save()
        for c in Candidate.objects.filter(sets=self):
            c.released = True
            c.save()


class Ballot(models.Model):
    """
    Ballots have manytomany fields pointing to a set of votes.
    Creating a new ballot takes a Set of candidates creates a
    vote that points to each and puts those votes on the ballot.
    """
    candidateset = models.ForeignKey(CandidateSet, null=True)
    name = models.CharField(max_length=255)
    voter = models.ForeignKey(User)
    creation_time = models.DateTimeField(default=timezone.now, null=True)
    completion_time = models.DateTimeField(null=True, blank=True)
    close_time = models.DateTimeField(null=True, blank=True)

    def closed(self):
        if self.close_time is None:
            return False
        else:
            return timezone.now() > self.close_time

    def import_candidateset(self, cset):
        self.candidateset = cset
        for candidate in Candidate.objects.filter(sets=cset):
            newvote = Vote(candidate=candidate, ballot=self)
            newvote.save()
        self.save()

    def remaining(self):
        return Vote.objects.filter(ballot=self, rating=None).count()

    def __unicode__(self):
        return self.name

VOTE_VALUES = {"up": 1, "down": -1, "meh": 0}

class Vote(models.Model):
    """A vote has a manytoone with a candidate. It has datetime, rating, comments, and user."""
    candidate = models.ForeignKey(Candidate)
    time = models.DateTimeField(null=True)
    rating = models.IntegerField(null=True)
    ballot = models.ForeignKey(Ballot, null=True)

    def rate(self, rating):
        """set vote rating.  if parent ballot is complete, set it's completion time"""
        self.rating = VOTE_VALUES[rating]
        self.time = timezone.now()
        self.save()
        if self.ballot.remaining() == 0:
            self.ballot.completion_time = timezone.now()
            self.ballot.save()

class Comment(models.Model):
    user = models.ForeignKey(User)
    candidate = models.ForeignKey(Candidate)
    time = models.DateTimeField()
    text = models.CharField(max_length=3000)

    def __unicode__(self):
        cutoff = 40
        if len(self.text) > cutoff:
            ellipse = "..."
        else:
            ellipse = ""
        return "{} - {}{} - {}".format(self.user.username, self.text[:cutoff], ellipse,  self.time)
