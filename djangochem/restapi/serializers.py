from django.contrib.auth.models import User, Group
from pgmols.models import *
from jobs.models import *
from rest_framework import serializers
from generic_relations.relations import GenericRelatedField

class UserSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = User
        fields = ('username', 'email', 'groups')


class GroupSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = Group
        fields = ('name',)



class JobConfigSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = JobConfig
        fields = ('name', 'configpath', 'parent_class_name')



class CompactJobConfigSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = JobConfig
        fields = ('name', 'url')



class MethodSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Method
        fields = ('name',
                  'description',
                  'details')


class MolSerializer(serializers.HyperlinkedModelSerializer):
    parentjob = serializers.HyperlinkedRelatedField(many=False, read_only=True, view_name='job-detail')
    childjobs = serializers.HyperlinkedRelatedField(many=True, read_only=True, view_name='job-detail')
    group = GroupSerializer(many=False, read_only=True)
    class Meta:
        model = Mol
        # lookup_field = "inchi_key"
        fields = ("id",
                  "inchikey",
                  "group",
                  "smiles",
                  "tags",
                  "mass",
                  "batches",
                  "createtime",
                  "parents",
                  "parentjob",
                  "childjobs",
                  "nicknames",
                  "details",
                  "url"
                  )


class CompactMolSerializer(serializers.HyperlinkedModelSerializer):
    group = GroupSerializer(many=False, read_only=True)
    class Meta:
        model = Mol
        fields = ("url",
                  "inchikey",
                  "group",
                  "smiles",
                  "tags",
                  "mass",
                  "batches",
                  "nicknames",
                  "details"
                  )



class CalcSerializer(serializers.HyperlinkedModelSerializer):
    parentjob = serializers.HyperlinkedRelatedField(many=False, read_only=True, view_name='job-detail')
    childjobs = serializers.HyperlinkedRelatedField(many=True, read_only=True, view_name='job-detail')
    mol = CompactMolSerializer(many=False, read_only=True)
    parents = serializers.HyperlinkedRelatedField(many=True, read_only=True, view_name='calc-detail')
    geoms = serializers.HyperlinkedRelatedField(many=True, read_only=True, view_name='geom-detail')
    method = MethodSerializer(many=False, read_only=True)

    class Meta:
        model = Calc
        fields = ('props',
                  'method',
                  'mol',
                  'parents',
                  'geoms',
                  'parentjob',
                  'childjobs')

class GeomSerializer(serializers.HyperlinkedModelSerializer):
    parentjob = serializers.HyperlinkedRelatedField(many=False, read_only=True, view_name='job-detail')
    childjobs = serializers.HyperlinkedRelatedField(many=True, read_only=True, view_name='job-detail')
    mol = CompactMolSerializer(many=False, read_only=True)
    parents = serializers.HyperlinkedRelatedField(many=True, read_only=True, view_name='geom-detail')
    method = MethodSerializer(many=False, read_only=True)

    class Meta:
        model = Geom
        fields = ('method',
                  'mol',
                  'parents',
                  'parentjob',
                  'childjobs',
                  'details',
                  'xyz'
                  )


class JobSerializer(serializers.HyperlinkedModelSerializer):
    parent = GenericRelatedField({
        Mol: CompactMolSerializer(),
        Geom: GeomSerializer(),
        Calc: CalcSerializer()
    })

    config = CompactJobConfigSerializer(many=False, read_only=True)
    class Meta:
        model = Job
        fields = ('config',
                  'group',
                  'priority',
                  'workbatch',
                  'status',
                  'createtime',
                  'details',
                  'completetime',
                  'duration',
                  'cores',
                  'architecture',
                  'repoversion',
                  'parent')
