from django.contrib.auth.models import User, Group
from django.db.models import F, Q, When, Value, Case
from rest_framework import viewsets
from pgmols.models import Mol, Geom, Calc, Method
from jobs.models import Job, JobConfig

from .serializers import UserSerializer, GroupSerializer
from .serializers import MolSerializer, JobSerializer, JobConfigSerializer, CalcSerializer, GeomSerializer, MethodSerializer
from rest_framework.response import Response
from rest_framework import permissions
from rest_framework.settings import api_settings
from rest_framework.pagination import LimitOffsetPagination

from rest_framework_csv import renderers as csv_renderers

from . import filters
from .svg_renderer import SVGRenderer

from six import text_type


class IsGroupMember(permissions.BasePermission):
    def has_object_permission(self, request, view, obj):
        user = request.user
        if user.is_staff:
            return True
        else:
            groups = user.groups.values_list('name', flat=True)
            return obj.released and obj.project in groups

    def has_permission(self, request, view):
        user = request.user
        if user.is_staff:
            return True
        else:
            return False


class HasMolPermissions(permissions.BasePermission):

    def has_permission(self, request, view):
        return True

    def has_object_permission(self, request, view, obj):
        user = request.user
        if user.is_staff:
            return True
        else:
            groups = user.groups.values_list('name', flat=True)
            return obj.mol.released and obj.mol.project in groups



class UserViewSet(viewsets.ModelViewSet):
    """
    API endpoint that allows users to be viewed or edited.
    """
    queryset = User.objects.all()
    serializer_class = UserSerializer


class GroupViewSet(viewsets.ModelViewSet):
    """
    API endpoint that allows groups to be viewed or edited.
    """
    queryset = Group.objects.all()
    serializer_class = GroupSerializer



class PaginatedCSVRenderer(csv_renderers.CSVRenderer):
    results_field = 'results'

    def render(self, data, media_type=None, renderer_context=None):
        if not isinstance(data, list):
            data = data.get(self.results_field, [])
        return super(PaginatedCSVRenderer, self).render(data, media_type, renderer_context)



class LargeResultsSetPagination(LimitOffsetPagination):
    default_limit = 100
    max_limit = 200000



# class MolCSVRenderer(PaginatedCSVRenderer):
#     """
#     Special Mol CSV renderer to handle 'properties' & 'sets' and output nice CSVs
#     """
#     results_field = 'results'

#     def flatten_dict(self, d):
#         flat_dict = {}
#         for key, item in d.items():
#             key = text_type(key)
#             if key == "properties":
#                 flat_item = self.flatten_properties(item)
#             elif key == "sets":
#                 flat_item = self.flatten_sets(item)
#             else:
#                 flat_item = self.flatten_item(item)
#             nested_item = self.nest_flat_item(flat_item, key)
#             flat_dict.update(nested_item)
#         return flat_dict

#     def flatten_properties(self, prop_list):
#         prop_dict = {}
#         for prop in prop_list:
#             prop_dict[prop["method"]+"."+prop["name"]] = prop["value"]
#         return prop_dict

#     def flatten_sets(self, sets):
#         set_dict = {}
#         for s in sets:
#             set_dict[s] = True
#         return set_dict


class MethodViewSet(viewsets.ModelViewSet):
    """
    API endpoint that allows groups to be viewed or edited.
    """
    queryset = Method.objects.all()
    serializer_class = MethodSerializer
    # filter_class = filters.MolFilter
    renderer_classes = api_settings.DEFAULT_RENDERER_CLASSES
    orderby_properties = ()
    orderby = ()

class MolViewSet(viewsets.ModelViewSet):
    """
    API endpoint that allows groups to be viewed or edited.
    """
    queryset = Mol.objects.all()
    serializer_class = MolSerializer
    # filter_class = filters.MolFilter
    renderer_classes = api_settings.DEFAULT_RENDERER_CLASSES + [SVGRenderer] # + [MolCSVRenderer]
    pagination_class = LargeResultsSetPagination
    permission_classes = [IsGroupMember]

    orderby_properties = ()
    orderby = ()

class GeomViewSet(viewsets.ModelViewSet):
    """
    API endpoint that allows groups to be viewed or edited.
    """
    queryset = Geom.objects.all()
    serializer_class = GeomSerializer
    # filter_class = filters.GeomFilter
    renderer_classes = api_settings.DEFAULT_RENDERER_CLASSES
    permission_classes = [IsGroupMember]

    orderby_properties = ()
    orderby = ()

class CalcViewSet(viewsets.ModelViewSet):
    """
    API endpoint that allows groups to be viewed or edited.
    """
    queryset = Calc.objects.all()
    serializer_class = CalcSerializer
    # filter_class = filters.CalcFilter
    renderer_classes = api_settings.DEFAULT_RENDERER_CLASSES
    pagination_class = LargeResultsSetPagination
    permission_classes = [IsGroupMember]

    orderby_properties = ()
    orderby = ()


class JobViewSet(viewsets.ModelViewSet):
    """
    API endpoint that allows groups to be viewed or edited.
    """
    queryset = Job.objects.all()
    serializer_class = JobSerializer


class JobConfigViewSet(viewsets.ModelViewSet):
    """
    API endpoint that allows groups to be viewed or edited.
    """
    queryset = JobConfig.objects.all()
    serializer_class = JobConfigSerializer


# class CommentViewSet(viewsets.ModelViewSet):
#     """
#     API endpoint that allows groups to be viewed or edited.
#     """
#     queryset = Comment.objects.all()
#     serializer_class = CommentSerializer
#     permission_classes = [HasMolPermissions]
#     filter_class = filters.CommentFilter

#     def get_queryset(self):
#         """
#         This view should return a list of all Mols viewable by the currently authenticated user
#         """
#         user = self.request.user
#         if user.is_staff:
#             return Comment.objects.all()
#         else:
#             groups = user.groups.values_list('name', flat=True)
#             return Comment.objects.filter(Mol__released=True, Mol__project__in=groups)


# class MethodViewSet(viewsets.ModelViewSet):
#     """
#     API endpoint that allows groups to be viewed or edited.
#     """
#     queryset = Method.objects.all()
#     serializer_class = MethodSerializer
#     filter_class = filters.MethodFilter

