from django.conf.urls import patterns, include, url
from django.contrib import admin
from django.conf import settings
from django.conf.urls.static import static
from django.contrib.auth.decorators import login_required, user_passes_test
from django.contrib.auth import views as auth_views


from molvote.views import BallotManageView, BallotResultsView

from django.views.generic.base import TemplateView

import registration
from registration.backends.default.views import ActivationView
from registration.backends.default.views import RegistrationView
from invitekey.views import InvitedRegistrationView, InviteView

admin.autodiscover()

urlpatterns = patterns('',
    # Examples:
    url(r'^$', 'molvote.views.home', name='home'),
    url(r'^samsung$', 'molvote.views.samsung_home', name='samsung_home'),
    url(r'^flow_batt$', 'molvote.views.flow_batt_home', name='flow_batt_home'),
    url(r'^laser$', 'molvote.views.laser_home', name='laser_home'),

    url(r'^table$', 'molvote.views.table', name='table'),

    url(r'^cards$', 'molvote.views.cards', name='cards'),
    url(r'^redox_cards$', 'molvote.views.redox', name='redox_cards'),
    url(r'^redox_bubble$', 'molvote.views.redox', name='redox_bubble'),
    url(r'^laser_cards$', 'molvote.views.laser_cards', name='laser_cards'),
    url(r'^laser_table$', 'molvote.views.laser_table', name='laser_table'),


    url(r'^card/(?P<inchi_key>[A-Z-]+)$', 'molvote.views.card', name='card'),
    url(r'^geometries/(?P<inchi_key>[A-Z-]+)$', 'molvote.views.geometries', name='geometries'),

    url(r'^detail/(?P<project_name>[a-zA-Z-]+)/(?P<inchi_key>[A-Z-]+)$', 'molvote.views.detail', name='detail'),

    url(r'^detail/(?P<candidate_id>[0-9]+)$', 'molvote.views.detail_by_id', name='detail_by_id'),


    url(r'^bubble$', 'molvote.views.bubble', name='bubble'),

    url(r'^redox_graph_json/(?P<redoxpair_id>[0-9]+)$', 'molvote.views.redox_graph_json', name='redox_graph_json'),
    url(r'^redox_graph/(?P<redoxpair_id>[0-9]+)$', 'molvote.views.redox_graph', name='redox_graph'),

    url(r'^redox_tree_json/(?P<redoxpair_id>[0-9]+)$', 'molvote.views.redox_tree_json', name='redox_tree_json'),
    url(r'^redox_tree/(?P<redoxpair_id>[0-9]+)$', 'molvote.views.redox_tree', name='redox_tree'),
    url(r'^super_tree/(?P<count>[0-9]+)$', 'molvote.views.redox_super_tree', name='redox_super_tree'),

    url(r'^statistics$', 'molvote.views.statistics', name='statistics'),

    url(r'^ballot_table/(?P<ballot_id>[0-9]+)$', 'molvote.views.ballot_table', name="ballot_table"),
    url(r'^ballots/(?P<ballot_id>[0-9]+)$', 'molvote.views.ballot', name="ballot"),
    url(r'^ballots/(?P<ballot_id>[0-9]+)/vote$', 'molvote.views.vote', name="vote"),

    url(r'^ballots/manage$', BallotManageView.as_view()),
    url(r'^ballots/results/(?P<set_id>[0-9]+)$', BallotResultsView.as_view(), name="results"),

    url(r'^chemworkers/', include('chemworkers.urls')),
    url(r'^api/', include('api.urls')),
    url(r'^rf/', include('rfapi.urls')),

    url(r'^admin/', include(admin.site.urls)),


    url(r'^accounts/invite', InviteView.as_view(), name="create_invitations"),

    url(r'^accounts/activate/complete/$',
       TemplateView.as_view(template_name='registration/activation_complete.html'),
       name='registration_activation_complete'),

    url(r'^accounts/activate/complete_preapproved/$',
       TemplateView.as_view(template_name='registration/activation_complete_preapproved.html'),
       name='registration_activation_complete_preapproved'),

    # Activation keys get matched by \w+ instead of the more specific
    # [a-fA-F0-9]{40} because a bad activation key should still get to the view;
    # that way it can return a sensible "invalid key" message instead of a
    # confusing 404.
    url(r'^accounts/activate/(?P<activation_key>\w+)/$',
       ActivationView.as_view(),
       name='registration_activate'),
    url(r'^accounts/register/$',
       InvitedRegistrationView.as_view(),
       name='registration_register'),
    url(r'^accounts/register/complete/$',
       TemplateView.as_view(template_name='registration/registration_complete.html'),
       name='registration_complete'),
    url(r'^accounts/register/closed/$',
       TemplateView.as_view(template_name='registration/registration_closed.html'),
       name='registration_disallowed'),
    url(r'^accounts/login/$',
       auth_views.login,
       {'template_name': 'registration/login.html'},
       name='login'),
    url(r'^accounts/logout/$',
       auth_views.logout,
       {'template_name': 'registration/logout.html'},
       name='logout'),
    url(r'^accounts/password/change/$',
       auth_views.password_change,
       name='password_change'),
    url(r'^accounts/password/change/done/$',
       auth_views.password_change_done,
       name='password_change_done'),
    url(r'^accounts/password/reset/$',
       auth_views.password_reset,
       name='password_reset'),
    url(r'^accounts/password/reset/confirm/(?P<uidb64>[0-9A-Za-z]+)-(?P<token>.+)/$',
       auth_views.password_reset_confirm,
       name='password_reset_confirm'),
    url(r'^accounts/password/reset/complete/$',
       auth_views.password_reset_complete,
       name='password_reset_complete'),
    url(r'^accounts/password/reset/done/$',
       auth_views.password_reset_done,
       name='password_reset_done'),


    url(r'^user/profile$', 'molvote.views.profile', name='profile'),

 )


if settings.DEBUG:
    import debug_toolbar
    urlpatterns += patterns('',
        (r'^__debug__/', include(debug_toolbar.urls)),
    )
    urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)