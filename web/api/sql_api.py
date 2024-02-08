import operator
import json
import csv
import StringIO

from django.db.models import Q
from django.core.serializers.json import DjangoJSONEncoder
from django.contrib.auth.models import User

from tastypie import fields
from tastypie.resources import ModelResource, ALL, ALL_WITH_RELATIONS
from tastypie.serializers import Serializer
from tastypie.cache import SimpleCache

from molvote.models import Candidate, CandidateSet, Vote, RedoxPair, Ballot
from security import ApiKeyPlusWebAuthentication, SQLAuthorization, VoteAuthorization
from blockassembler import blockdraw
from blocks import block

from django.conf import settings

IMG_CACHE = settings.IMG_CACHE_DIR

ROUND_TO = {"weight": 2,
            "sascore": 1,
            "splitting": 3,
            "strength": 3,
            "absorption": 2,
            "rate": 2,
            "homo": 2,
            "lumo": 2,
            "water_solvation_energy": 2
            }


class CandidateSerializer(Serializer):
    formats = ['json', 'jsonp', 'xml', 'yaml', 'html', 'plist', 'smi', 'svg', 'csv']
    content_types = {
        'json': 'application/json',
        'jsonp': 'text/javascript',
        'xml': 'application/xml',
        'yaml': 'text/yaml',
        'html': 'text/html',
        'plist': 'application/x-plist',
        'smi': 'text/plain',
        'svg': 'image/svg+xml',
        'csv': 'text/csv'
    }

    def to_smi(self, data, options=None):
        options = options or {}
        data = self.to_simple(data, options)
        if "objects" in data:
            return "\n".join([i["smiles"] for i in data["objects"]]) + "\n"
        else:
            return data["smiles"] + "\n"

    def to_html(self, data, options=None):
        options = options or {}
        data = self.to_simple(data, options)

        if "objects" in data and len(data["objects"]) > 0:
            out = "\n".join(["<a href='{0[resource_uri]}'><img width=200px height=200px src='{0[resource_uri]}?format=svg'></a>".format(i) for i in data["objects"]])
            return out

    def to_json(self, data, options=None):
        options = options or {}

        data = self.to_simple(data, options)
        if "objects" in data:
            for d in data["objects"]:
                for k in d:
                    if k in ROUND_TO and k in d and d[k]:
                        d[k] = round(d[k], ROUND_TO[k])
        else:
            for k in data:
                if k in ROUND_TO and k in data and data[k]:
                    data[k] = round(data[k], ROUND_TO[k])
        return json.dumps(data, cls=DjangoJSONEncoder, sort_keys=True)


    def to_csv(self, data, options=None):
        options = options or {}
        data = self.to_simple(data, options)
        if "objects" in data:
            output = StringIO.StringIO()
            writer = csv.DictWriter(output, data["objects"][0].keys())
            writer.writeheader()
            writer.writerows(data["objects"])
            return output.getvalue()


    def to_svg(self, data, options=None):
        options = options or {}
        data = self.to_simple(data, options)
        if "objects" in data:
            data = data["objects"][0]
        try:
            filename = blockdraw.write_svg(block.Block(data["smiles"]), path=IMG_CACHE, create_groups=data["project"])
        except:
            filename = blockdraw.write_svg(block.Block(data["smiles"]), path=IMG_CACHE, molconvert=True)

        return open(filename, 'r').read()


class CandidateResource(ModelResource):
    class Meta:
        authentication = ApiKeyPlusWebAuthentication()
        authorization = SQLAuthorization()
        serializer = CandidateSerializer()
        limit = 100
        max_limit = 1000
        queryset = Candidate.objects.all()
        allowed_methods = ['get']
        ordering = ["inchi_key",
                    "nicknames",
                    "absorption",
                    "splitting",
                    "strength",
                    "rate",
                    "weight",
                    "sascore",
                    "donor",
                    "acceptor",
                    "bridge1",
                    "bridge2",
                    "homo",
                    "lumo",
                    "vote",
                    "calc_time"]
        filtering = {"project": ALL,
                    "vote": ALL_WITH_RELATIONS,
                    "sets": ALL_WITH_RELATIONS,
                    "smiles": ALL,
                    "inchi_key": ALL,
                    "smiles": ALL,
                    "absorption": ALL,
                    "splitting": ALL,
                    "strength": ALL,
                    "rate": ALL,
                    "homo": ALL,
                    "lumo": ALL,
                    "nicknames": ALL,
                    "calc_time": ALL,
                    "sascore": ALL,
                    "weight": ALL,
                    "water_solvation_energy": ALL,
                    "total_energy": ALL,
                    "cas": ALL,
                    "donor":ALL_WITH_RELATIONS,
                    "acceptor":ALL_WITH_RELATIONS,
                    "bridge1":ALL_WITH_RELATIONS,
                    "bridge2":ALL_WITH_RELATIONS,}

    sets = fields.ToManyField("api.sql_api.CandidateSetResource", "sets", null=True, full=False)
    donor = fields.ToOneField("api.sql_api.CandidateResource", "donor", null=True, full=False)
    acceptor = fields.ToOneField("api.sql_api.CandidateResource", "acceptor", null=True, full=False)
    bridge1 = fields.ToOneField("api.sql_api.CandidateResource", "bridge1", null=True, full=False)
    bridge2 = fields.ToOneField("api.sql_api.CandidateResource", "bridge2", null=True, full=False)


    def build_filters(self, filters=None):
        if filters is None:
            filters = {}

        special_filters = {}
        for k in filters.keys():
            if k.startswith("not__"):
                special_filters.update({k: filters.pop(k)})
            if k.endswith("__all"):
                special_filters.update({k: filters.pop(k)})

        orm_filters = super(CandidateResource, self).build_filters(filters)
        orm_filters.update(special_filters)
        return orm_filters

    def apply_filters(self, request, applicable_filters):
        if applicable_filters is None:
            applicable_filters = {}
        negated = {}
        unions = []
        for k in applicable_filters.keys():
            if k.startswith("not__"):
                negated[k[5:]] = applicable_filters.pop(k)[0].split(',')
            if k.endswith("__all"):
                unions.append((k[:-5], applicable_filters.pop(k)[0].split(',')))

        semi_filtered = super(CandidateResource, self).apply_filters(request, applicable_filters)

        for union in unions:
            print union[1]
            for val in union[1]:
                union_query = {union[0]:val}
                semi_filtered = semi_filtered.filter(**union_query)
        if negated:
            semi_filtered = semi_filtered.exclude(**negated)
        return semi_filtered


class RedoxPairResource(ModelResource):
    class Meta:
        authentication = ApiKeyPlusWebAuthentication()
        authorization = SQLAuthorization()
        limit = 100
        max_limit = 5000
        queryset = RedoxPair.objects.all()
        allowed_methods = ['get']
        ordering = ["redox_potential", "log_hyd_constant", "reduced", "substituent_count"]
        filtering = {"reduced": ALL_WITH_RELATIONS,
                    "oxidized": ALL_WITH_RELATIONS,
                    "hydrated_oxidized": ALL_WITH_RELATIONS,
                    "log_hyd_constant": ALL,
                    "log_hyd_constant_error": ALL,
                    "redox_potential": ALL,
                    "redox_potential_error": ALL,
                    "substituent_count": ALL,
                    "is_minimum": ALL,
                    "michael_hydration_energy": ALL,
                     }
        cache = SimpleCache(timeout=3600)

    reduced = fields.ToOneField("api.sql_api.CandidateResource", "reduced", full=True)
    oxidized = fields.ToOneField("api.sql_api.CandidateResource", "oxidized", full=True)
    hydrated_oxidized = fields.ToOneField("api.sql_api.CandidateResource", "hydrated_oxidized", null=True, full=True)



    def apply_sorting(self, obj_list, options=None):
        """
        Allows for the sorting of objects being returned.
        This needs to be implemented at the user level.
        ``ModelResource`` includes a full working version specific to Django's
        ``Models``.
        """
        if "raw_order_by" in options:
            return obj_list.order_by(options["raw_order_by"])
        return obj_list

    def build_filters(self, filters=None):
        if filters is None:
            filters = {}

        special_filters = {}
        for k in filters.keys():
            if k.startswith("not__"):
                special_filters.update({k: filters.pop(k)})
            if k.endswith("__all"):
                special_filters.update({k: filters.pop(k)})

        orm_filters = super(RedoxPairResource, self).build_filters(filters)
        orm_filters.update(special_filters)
        return orm_filters

    def apply_filters(self, request, applicable_filters):
        if applicable_filters is None:
            applicable_filters = {}
        negated = {}
        unions = []
        for k in applicable_filters.keys():
            if k.startswith("not__"):
                negated[k[5:]] = applicable_filters.pop(k)[0].split(',')
            if k.endswith("__all"):
                unions.append((k[:-5], applicable_filters.pop(k)[0].split(',')))

        semi_filtered = super(RedoxPairResource, self).apply_filters(request, applicable_filters)

        for union in unions:
            print union[1]
            for val in union[1]:
                union_query = {union[0]:val}
                semi_filtered = semi_filtered.filter(**union_query)
        if negated:
            semi_filtered = semi_filtered.exclude(**negated)

        if "all" not in applicable_filters or applicable_filters["all"] == False:
            semi_filtered = semi_filtered.filter(is_minimum=True)

        return semi_filtered


class MiniRedoxPairResource(RedoxPairResource):
    class Meta:
        authentication = ApiKeyPlusWebAuthentication()
        authorization = SQLAuthorization()
        limit = 10000
        max_limit = 50000
        queryset = RedoxPair.objects.all()
        allowed_methods = ['get']
        ordering = ["redox_potential", "log_hyd_constant", "reduced", "substituent_count"]
        filtering = {"reduced": ALL_WITH_RELATIONS,
                    "oxidized":ALL_WITH_RELATIONS,
                    "hydrated_oxidized":ALL_WITH_RELATIONS,
                    "log_hyd_constant": ALL,
                    "redox_potential": ALL,
                    "substituent_count": ALL,
                    "is_minimum": ALL,
                    }
        fields = ["log_hyd_constant", "redox_potential", "water_solvation_energy", "substituent_count"]
    reduced = fields.ToOneField("api.sql_api.RedoxCandidateResource", "reduced", full=True, use_in='detail')
    oxidized = fields.ToOneField("api.sql_api.RedoxCandidateResource", "oxidized", full=True, use_in='detail')
    hydrated_oxidized = fields.ToOneField("api.sql_api.RedoxCandidateResource", "hydrated_oxidized", null=True, full=True, use_in='detail')

class RedoxCandidateResource(CandidateResource):
    class Meta:
        authentication = ApiKeyPlusWebAuthentication()
        authorization = SQLAuthorization()
        queryset = Candidate.objects.all()
        allowed_methods = ['get']
        ordering = ["water_solvation_energy"]
        fields = ["inchi_key","id","water_solvation_energy"]
        filtering = {"inchi_key": ALL,
                     "water_solvation_energy": ALL,
                     "sets": ALL_WITH_RELATIONS}
        limit = 500
        max_limit = 5000
        serializer = CandidateSerializer()

    sets = fields.ToManyField("api.sql_api.CandidateSetResource", "sets", full=False, use_in='detail')


class MiniCandidateResource(CandidateResource):
    class Meta:
        authentication = ApiKeyPlusWebAuthentication()
        authorization = SQLAuthorization()
        queryset = Candidate.objects.all()
        allowed_methods = ['get']
        ordering = ["absorption", "splitting", "strength", "rate", "weight", "sascore"]
        fields = ["splitting", "strength", "rate", "absorption", "weight", "sascore", "inchi_key","id","water_solvation_energy"]
        filtering = {"sets": ALL_WITH_RELATIONS,
                    "smiles": ALL,
                    "inchi_key": ALL,
                    "absorption": ALL,
                    "splitting": ALL,
                    "strength": ALL,
                    "rate": ALL,
                    "homo": ALL,
                    "lumo": ALL,
                    "nicknames": ALL,
                    "calc_time": ALL,
                    "sascore": ALL,
                    "weight": ALL,
                    "water_solvation_energy": ALL}
        limit = 5000
        max_limit = 100000
        serializer = CandidateSerializer()

    sets = fields.ToManyField("api.sql_api.CandidateSetResource", "sets", full=False, use_in='detail')
    calc_time = fields.DateTimeField(use_in='detail', attribute='calc_time')


class CandidateSetResource(ModelResource):
    class Meta:
        authentication = ApiKeyPlusWebAuthentication()
        authorization = SQLAuthorization()
        queryset = CandidateSet.objects.all()
        allowed_methods = ['get']
        filtering = {"name": ALL}


class BallotResource(ModelResource):
    class Meta:
        authentication = ApiKeyPlusWebAuthentication()
        queryset = Ballot.objects.all()
        allowed_methods = ['get']
        filtering = {"candidateset": ALL_WITH_RELATIONS,
                     "voter": ALL_WITH_RELATIONS}

    candidateset = fields.ToOneField("api.sql_api.CandidateSetResource", "candidateset", full=False, use_in='detail')
    voter = fields.ToOneField("api.sql_api.UserResource", "voter", full=False)


class VoteResource(ModelResource):
    class Meta:
        authentication = ApiKeyPlusWebAuthentication()
        authorization = VoteAuthorization()
        queryset = Vote.objects.all()
        allowed_methods = ['get']
        ordering = ["rating"]
        filtering = {"ballot": ALL_WITH_RELATIONS,
                     "candidate": ALL}
        limit = 1000
        max_limit = 10000

    ballot = fields.ToOneField("api.sql_api.BallotResource", "ballot", full=False)
    candidate = fields.ToOneField("api.sql_api.CandidateResource", "candidate", full=False)


class UserResource(ModelResource):
    class Meta:
        queryset = User.objects.all()
        authentication = ApiKeyPlusWebAuthentication()
        allowed_methods = ['get']
        fields = ["username"]
        filtering = {"username": ALL,
                     "id": ALL}
