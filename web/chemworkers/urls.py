from django.conf.urls import patterns, url

urlpatterns = patterns('',
    url(r'^latest/?$',  'chemworkers.views.latest', name="chemworkers_latest_all"),
    url(r'^latest/(?P<worker_names>[,+_A-Za-z0-9-]+)$',  'chemworkers.views.latest', name="chemworkers_latest")
)
