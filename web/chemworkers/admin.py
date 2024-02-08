from django.contrib import admin
from django import forms
from django.conf import settings

from models import WorkLog, WorkLocation, ChemWorker
# Register your models here.

from workers import configs

import fnmatch
import os
from operator import itemgetter

config_base = getattr(settings, "CHEM_WORKERS_CONFIG_BASE")

class ChemWorkerForm(forms.ModelForm):
    def __init__(self, *args, **kwargs):
        forms.ModelForm.__init__(self, *args, **kwargs)

        worker_choices = []
        base_len = len(config_base)
        for root, dirnames, filenames in os.walk(config_base):
            for filename in fnmatch.filter(filenames, '*.json'):
                if not root.endswith("dbs") and not root.endswith("test"):
                    worker_choices.append(os.path.join(root[base_len:], filename).strip("/"))


        #labels = ["/".join(path.split("/")[-2:]) for path in worker_choices]
        pairs = zip(worker_choices, worker_choices)
        pairs.sort()  # key=itemgetter(1)
        self.fields["config_path"] = forms.ChoiceField(choices=pairs)

    class Meta:
        model = ChemWorker
        fields = ["name",
                  "label",
                  "config_path",
                  "project",
                  "work_location",
                  "batch_size",
                  "request_count_min",
                  "request_count_max",
                  
                  "min_priority",
                  "max_priority",
                  "purge_mode",
                  "ignore_lock",

                  ]


class ChemWorkerAdmin(admin.ModelAdmin):
    form = ChemWorkerForm

admin.site.register(WorkLog)
admin.site.register(WorkLocation)
admin.site.register(ChemWorker, ChemWorkerAdmin)
