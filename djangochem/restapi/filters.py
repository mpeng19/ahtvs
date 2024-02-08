from pgmols.models import *

# import rest_framework_filters as filters
import django_filters

from django.forms import Field
from django_filters.filters import Filter, DateTimeFilter
from django.db.models import AutoField
from django.db.models.constants import LOOKUP_SEP

from rest_framework import filters


STRING_LOOKUPS = ['exact', 'contains', 'in', 'startswith']
NUMBER_LOOKUPS = ['exact', 'lt', 'lte', 'gt', 'gte']


class IsOwnerFilterBackend(filters.BaseFilterBackend):
    """
    Filter that only allows users to see their own objects.
    """
    def filter_queryset(self, request, queryset, view):
        return queryset.filter(owner=request.user)


# class CandidateSetFilter(django_filters.FilterSet):
#     name = django_filters.CharFilter(name="name", lookup_type="contains")

#     class Meta:
#         model = CandidateSet
#         fields = ["name"]

PROPERTY_ATT_NAME = "properties"
NAME_ATT = "name"
VALUE_ATT = "value"
NAME_KEY = PROPERTY_ATT_NAME+LOOKUP_SEP+NAME_ATT


class CommaField(Field):
    def clean(self, value):  # convert 'this,that,other' to {'this', 'that', 'other'}
        return set(v for v in value.split(',')) if value else ()


class CommaInFilter(django_filters.CharFilter):
    field_class = CommaField


class CommaAllFilter(django_filters.MultipleChoiceFilter):
    field_class = CommaField


class CaseInsensitiveBooleanFilter(django_filters.Filter):
    def filter(self, qs, value):
        if value is not None:
            lc_value = value.lower()
            if lc_value == "true":
                value = True
            elif lc_value == "false":
                value = False
            return qs.filter(**{self.name: value})
        return qs


# class CandidateFilter(django_filters.FilterSet):
#     sets = django_filters.CharFilter(name='sets__name', lookup_type="exact")
#     sets__in = CommaInFilter(name='sets__name', lookup_type="in")
#     sets__not = django_filters.CharFilter(name='sets__name', lookup_type="exact", exclude=True)
#     sets__contains = django_filters.CharFilter(name='sets__name', lookup_type="contains")
#     sets__all = CommaAllFilter(name='sets__name', conjoined=True)
#     released = CaseInsensitiveBooleanFilter(name='released', lookup_type="exact")

#     votes = django_filters.NumberFilter(name='vote__rating', lookup_type="exact")
#     votes__gte = django_filters.NumberFilter(name='vote__rating', lookup_type="gte")
#     votes__lte = django_filters.NumberFilter(name='vote__rating', lookup_type="lte")
#     votes__lt = django_filters.NumberFilter(name='vote__rating', lookup_type="lt")
#     votes__gt = django_filters.NumberFilter(name='vote__rating', lookup_type="gt")
#     prop = django_filters.MethodFilter(action='filter_property', lookup_type=None)

#     inchi_key = django_filters.CharFilter(name='inchi_key', lookup_type="exact")
#     inchi_key__contains = django_filters.CharFilter(name='inchi_key', lookup_type="contains")
#     inchi_key__in = CommaInFilter(name='inchi_key', lookup_type="in")

#     donor = django_filters.ModelChoiceFilter(queryset=Candidate.objects.all())
#     acceptor = django_filters.ModelChoiceFilter(queryset=Candidate.objects.all())

#     class Meta:
#         model = Candidate
#         fields = {'sets': ["exact", "in"],
#                   'released': ['exact'],
#                   'project': STRING_LOOKUPS,
#                   'votes': NUMBER_LOOKUPS,
#                   'smiles': STRING_LOOKUPS,
#                   'nicknames': STRING_LOOKUPS,
#                   'calc_time': NUMBER_LOOKUPS,
#                   'sascore': NUMBER_LOOKUPS,
#                   'weight': NUMBER_LOOKUPS,
#                   'prop': None}

#         # order_by = ['weight', 'sascore']
#     def filter_property(self, queryset, value):
#         method = self.data.get("method", "b3lyp_tddft_631gs_rpa_s0_geom")
#         for key, value in self.data.items():
#             if key.startswith('prop-'):
#                 name_lookup = key[5:].split(LOOKUP_SEP)
#                 name = name_lookup[0]
#                 if len(name_lookup) > 1:
#                     valuekey = PROPERTY_ATT_NAME+LOOKUP_SEP+VALUE_ATT+LOOKUP_SEP+name_lookup[1]
#                 else:
#                     valuekey = PROPERTY_ATT_NAME+LOOKUP_SEP+VALUE_ATT
#                 filterdict = {NAME_KEY: name, valuekey: value}
#                 if method:
#                     filterdict[PROPERTY_ATT_NAME+LOOKUP_SEP+"method__name"] = method
#                 queryset = queryset.filter(**filterdict)
#         return queryset


# class RedoxPairFilter(django_filters.FilterSet):
#     sets = django_filters.CharFilter(name="reduced__sets__name", lookup_type=["exact", "contains", "in"])

#     class Meta:
#         model = RedoxPair
#         fields = ["sets"]


# class IDField(Field):
#     def clean(self, value):  # convert '1,2,3' to {1, 2, 3}
#         return set(int(v) for v in value.split(',') if v.isnumeric()) if value else ()


# class IDFilter(Filter):
#     field_class = IDField


# class CommentFilter(django_filters.FilterSet):
#     candidate = IDFilter(name="candidate", lookup_type="in")
#     order_by_field = "ordering"

#     class Meta:
#         model = Comment
#         fields = {'id': ('in', 'exact'),
#                   'candidate': ('in', 'exact'),
#                   'time': ('gte', 'lte', 'exact')
#                   }
#         order_by = True


# class MethodFilter(django_filters.FilterSet):
#     name = django_filters.CharFilter(name="name", lookup_type=["exact", "contains", "in"])

#     class Meta:
#         model = Method
#         fields = ["name"]


# class PropertyFilter(django_filters.FilterSet):
#     candidate = django_filters.ModelChoiceFilter(queryset=Candidate.objects.all())

#     class Meta:
#         model = Property
#         fields = ["candidate"]
