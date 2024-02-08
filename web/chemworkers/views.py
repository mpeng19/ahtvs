import datetime
import collections

from django.shortcuts import render
from django.utils import timezone

from utils.loginhelpers import only_staff_allowed
from .models import WorkLog, ChemWorker
# Create your views here.


@only_staff_allowed
def latest(request, worker_names=None):
    now = timezone.now()
    if worker_names:
        logs = WorkLog.objects.filter(timestamp__gt=now-datetime.timedelta(days=30), chem_worker__name__in=worker_names.split(","))
        title = logs[0].chem_worker.config_path.split("/")[0]
    else:
        title = "All"
        logs = WorkLog.objects.filter(timestamp__gt=now-datetime.timedelta(days=30))
    # round now to next hour
    now = now - datetime.timedelta(minutes=now.minute,
                         seconds=now.second,
                         microseconds=now.microsecond)
    now = now + datetime.timedelta(hours=1)

    lastday = [l for l in logs if now - l.timestamp < datetime.timedelta(1)]


    hours = [[ (now - datetime.timedelta(hours=i+1)), 0, 0, 0, 0] for i in range(24)]
    for l in lastday:
        hours_ago = int((now - l.timestamp).seconds / 3600.)
        hour = hours[hours_ago]
        hour[1] += l.jobs_requested
        hour[2] += l.jobs_completed
        hour[3] += l.jobs_failed
        hour[4] = hour[2] + hour[3]

    now = now - datetime.timedelta(hours=now.hour,
                         minutes=now.minute,
                         seconds=now.second,
                         microseconds=now.microsecond)
    now = now + datetime.timedelta(days=1)
    last30 = [l for l in logs if now - l.timestamp < datetime.timedelta(30)]

    days = [[(now - datetime.timedelta(days=i)), 0, 0, 0, 0] for i in range(30)]
    for l in last30:
        days_ago = int((now - l.timestamp).days)
        day = days[days_ago]
        day[1] += l.jobs_requested
        day[2] += l.jobs_completed
        day[3] += l.jobs_failed
        day[4] = day[2] + day[3]


    max_hour_val = max([max(h[1:]) for h in hours])
    max_day_val = max([max(h[1:]) for h in days])

    worker_groups = collections.OrderedDict()
    for w in ChemWorker.objects.all().order_by("config_path"):
        label = w.config_path.split("/")[0]
        group = worker_groups.get(label, [])
        group.append(w.name)
        worker_groups[label] = group

    context = {'title':title,
                'hours': hours,
                'days': days,
                'max_hour_val': max_hour_val,
                'max_day_val': max_day_val,
                'worker_groups': worker_groups}
    return render(request, 'chemworkers/latest.html', context)

