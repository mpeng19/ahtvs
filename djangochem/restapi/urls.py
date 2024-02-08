from django.conf.urls import url, include
from rest_framework import routers
from . import views

router = routers.DefaultRouter()
router.register(r'users', views.UserViewSet)
router.register(r'groups', views.GroupViewSet)
router.register(r'mols', views.MolViewSet)
router.register(r'calcs', views.CalcViewSet)
router.register(r'geoms', views.GeomViewSet)
router.register(r'methods', views.MethodViewSet)
router.register(r'jobs', views.JobViewSet)
router.register(r'jobconfigs', views.JobConfigViewSet)
# router.register(r'comments', views.CommentViewSet)
# router.register(r'candidatesets', views.CandidateSetViewSet)
# router.register(r'redoxpairs', views.RedoxPairViewSet)
# router.register(r'properties', views.PropertyViewSet)
# router.register(r'methods', views.MethodViewSet)

# Wire up our API using automatic URL routing.
# Additionally, we include login URLs for the browsable API.
urlpatterns = [
    url(r'^', include(router.urls)),
    url(r'^api-auth/', include('rest_framework.urls', namespace='rest_framework'))
]