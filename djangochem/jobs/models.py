import os
import uuid

from django.db import models
from django.contrib.postgres.fields import JSONField
from django.contrib.auth.models import User, Group
from django.utils import timezone

from django.contrib.contenttypes.fields import GenericForeignKey
from django.contrib.contenttypes.models import ContentType

PARENT_CLASS_CHOICES = (('Mol', 'Mol'),
                        ('Geom', 'Geom'),
                        ('Calc', 'Calc'))


class JobConfig(models.Model):
    name = models.CharField(max_length=128)
    configpath = models.CharField(max_length=256)
    parent_class_name = models.CharField(max_length=4, choices=PARENT_CLASS_CHOICES)

    @property
    def dir(self):
        return os.path.dirname(self.configpath)

    def __repr__(self):
        return '<JobConfig ' + self.name + '>'

    def __str__(self):
        return self.name

STATUS_CHOICES = (('claimed', 'Claimed'),
                  ('done', 'Done'),
                  ('error', 'Error'),
                  ('', 'None'))


class Job(models.Model):
    uuid = models.UUIDField(unique=True, default=uuid.uuid4)
    config = models.ForeignKey(JobConfig, on_delete=models.CASCADE, null=True) #TODO:check on delete
    group = models.ForeignKey(Group, on_delete=models.CASCADE, null=True) #TODO:check on delete
    priority = models.IntegerField(default=0, db_index=True)
    workbatch = models.ForeignKey('WorkBatch', on_delete=models.CASCADE, null=True) #TODO:check on delete
    status = models.CharField(max_length=24, blank=True, default='', choices=STATUS_CHOICES, db_index=True)
    createtime = models.DateTimeField(default=timezone.now)
    parentct = models.ForeignKey(ContentType,
                                 on_delete=models.CASCADE,
                                 related_name='childjob',
                                 blank=True,
                                 null=True)
    parentid = models.PositiveIntegerField(blank=True, null=True)
    parent = GenericForeignKey('parentct', 'parentid')
    details = JSONField(null=True)
    completetime = models.DateTimeField(null=True)
    duration = models.FloatField(null=True)
    cores = models.IntegerField(null=True)
    architecture = models.CharField(max_length=24, null=True)
    repoversion = models.CharField(max_length=24, null=True)

    def __str__(self):
        return str(self.id) + ' : ' + str(self.createtime)


class WorkBatch(models.Model):
    user = models.ForeignKey(User, on_delete=models.CASCADE, null=True) #TODO:check on delete
    createtime = models.DateTimeField(auto_now=True)
    comments = models.TextField()
    name = models.CharField(max_length=128)
    details = JSONField(null=True)

    def __str__(self):
        return str(self.id) + ' : ' + str(self.createtime)