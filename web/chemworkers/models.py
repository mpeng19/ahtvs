import os

from django.db import models
from django.utils import timezone
# Create your models here.


class WorkLocation(models.Model):
    label = models.CharField(max_length=255)
    request_path = models.CharField(max_length=512)
    complete_path = models.CharField(max_length=512)
    archive_path = models.CharField(max_length=512)
    error_path = models.CharField(max_length=512)

    def __unicode__(self):
        return self.label


class ChemWorker(models.Model):
    name = models.CharField(max_length=255, unique=True)
    label = models.CharField(max_length=256)
    config_path = models.CharField(max_length=512)  # path relative to root of samsung git repo
    ignore_lock = models.BooleanField(default=False)
    request_count_min = models.IntegerField(default=0)
    request_count_max = models.IntegerField(default=10)
    purge_mode = models.BooleanField(default=False)
    min_priority = models.IntegerField(default=0)
    max_priority = models.IntegerField(default=10)
    batch_size = models.IntegerField(default=1)
    project = models.CharField(max_length=256)
    work_location = models.ForeignKey(WorkLocation)

    def __unicode__(self):
        args = (self.label,
                self.name,
                os.path.basename(os.path.dirname(self.config_path)),
                os.path.basename(self.config_path), self.work_location.request_path)
        return "{} [{}] ({}/{} -> {})".format(*args)


class WorkLog(models.Model):
    timestamp = models.DateTimeField(default=timezone.now)
    chem_worker = models.ForeignKey(ChemWorker)
    jobs_requested = models.IntegerField(default=0)
    jobs_completed = models.IntegerField(default=0)
    jobs_failed = models.IntegerField(default=0)





# def create_task(w, interval):
#     new_task = PeriodicTask(name=w.name,
#                           task=u'chemworkers.do_one_round',
#                           enabled=False,
#                           interval=interval,
#                           args=u'["{}"]'.format(w.name))
#     new_task.save()


# def create_from_existing(w, newname):
#     new_cw = ChemWorker(name=w.name[:-1] + newname[0],
#              label=" ".join(w.label.split()[:-1])+" "+newname+"-od",
#              config_path=w.config_path,
#              project=w.project,
#              work_location= WorkLocation.objects.get(label__startswith=newname))
#     return new_cw
