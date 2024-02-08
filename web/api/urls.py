from django.conf.urls import patterns, include
from tastypie.api import Api

from api import MoleculeResource, ReactiveMolResource, MetaResource, CalculationResource, MolecularLinkageResource
from sql_api import CandidateResource, CandidateSetResource, MiniCandidateResource
from sql_api import VoteResource, RedoxPairResource, BallotResource, MiniRedoxPairResource

api = Api(api_name="latest")
api.register(MoleculeResource())
api.register(ReactiveMolResource())
api.register(MetaResource())
api.register(CalculationResource())
api.register(MolecularLinkageResource())

api.register(VoteResource())
api.register(BallotResource())
api.register(CandidateSetResource())
api.register(CandidateResource())
api.register(MiniCandidateResource())
api.register(RedoxPairResource())
api.register(MiniRedoxPairResource())


urlpatterns = patterns('',
    (r'^', include(api.urls))
)
