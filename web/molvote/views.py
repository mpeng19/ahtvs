import json
import re
import functools

from django.utils import timezone
from django.shortcuts import render
from django.http import HttpResponse
from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.shortcuts import get_object_or_404
from django.core.exceptions import PermissionDenied
from django.views.generic import View
from django.db.models import Max, Min

from utils.loginhelpers import StaffOnlyMixin, only_staff_allowed
from utils.loginhelpers import LoginRequiredMixin

from .models import Candidate, CandidateSet, Ballot, Vote, RedoxPair, Comment, Property
from .forms import NewBallotForm

from aag_python.molecular_storage import molecular_data_models as mdm


# ########  HELPER FUNCTIONS #######################


def count_votes(ballot):
    remaining = Vote.objects.filter(ballot=ballot, rating=None).count()
    meh = Vote.objects.filter(ballot=ballot, rating=0).count()
    up = Vote.objects.filter(ballot=ballot, rating__gt=0).count()
    down = Vote.objects.filter(ballot=ballot, rating__lt=0).count()
    total = remaining + meh + up + down

    tallies = {"total": total,
               "up": up,
               "down": down,
               "meh": meh,
               "remaining": remaining}
    return tallies


def get_sets_for_user(user, limit_groups=[]):
    if user.is_staff:
        if limit_groups:
            query = CandidateSet.objects.filter(project__in=limit_groups)
        else:
            query = CandidateSet.objects.all()
        set_names = [s.name for s in query.only("name")]
    else:
        groups = user.groups.values_list('name', flat=True)
        if limit_groups:
            groups = list(set(groups) & set(limit_groups))
        query = CandidateSet.objects.filter(released=True, project__in=groups)
        set_names = [s.name for s in query.only("name")]

    # def mycmp(a, b):
    #     res = cmp(len(a), len(b))
    #     if res != 0:
    #         return res
    #     else:
    #         return cmp(a, b)
    set_names.sort()
    return set_names

# ################# VIEWS #########################


@login_required
def profile(request):
    return render(request, template_name='profile.html')


@login_required
def generic_home(request, project_name):
    user_groups = request.user.groups.values_list('name', flat=True)
    if project_name not in user_groups:
        raise PermissionDenied()

    latest_sets = CandidateSet.objects.filter(released=True,
                                              announced=True,
                                              project=project_name)
    if request.user.is_staff:
        latest_sets = latest_sets.order_by('-creation_time')[:5]
        user_sets = CandidateSet.objects.filter(announced=False, project=project_name)
        user_sets = user_sets.order_by('-creation_time')[:5]
    else:

        latest_sets = latest_sets.order_by('-creation_time')[:5]
        user_sets = CandidateSet.objects.filter(creator=request.user,
                                                project__in=user_groups)
        user_sets = user_sets.order_by('-creation_time')[:5]

    context = {"latest_sets": latest_sets,
               "user_sets": user_sets,
               }
    return render(request, 'home_{}.html'.format(project_name), context)



@login_required
def flow_batt_home(request):
    return generic_home(request, "flow_batt")


@login_required
def laser_home(request):
    return generic_home(request, "laser")


@login_required
def samsung_home(request):
    user_groups = request.user.groups.values_list('name', flat=True)
    if "samsung" not in user_groups:
        raise PermissionDenied()

    unfinished_ballots = Ballot.objects.filter(voter=request.user, completion_time__isnull=True)

    latest_sets = CandidateSet.objects.filter(released=True,
                                              announced=True,
                                              project="samsung")
    if request.user.is_staff:
        latest_sets = latest_sets.order_by('-creation_time')[:5]
        user_sets = CandidateSet.objects.filter(announced=False).order_by('-creation_time')[:5]
        ballots = Ballot.objects.all()
        sets_with_ballots = set(b.candidateset for b in ballots if b.candidateset_id is not None)
    else:

        latest_sets = latest_sets.order_by('-creation_time')[:5]
        user_sets = CandidateSet.objects.filter(creator=request.user,
                                                project__in=user_groups)
        user_sets = user_sets.order_by('-creation_time')[:5]
        ballots = Ballot.objects.filter(close_time__lt=timezone.now(),
                                        candidateset__released=True,
                                        candidateset__project__in=user_groups)
        sets_with_ballots = set(b.candidateset for b in ballots if b.candidateset_id is not None)

    context = {"unfinished_ballots": unfinished_ballots,
               "latest_sets": latest_sets,
               "user_sets": user_sets,
               "sets_with_ballots": sets_with_ballots
               }
    return render(request, 'home.html', context)


@login_required
def card(request, inchi_key):
    can = get_object_or_404(Candidate, inchi_key=inchi_key)
    context = {"candidate": can}
    return render(request, 'card.html', context)


@login_required
def home(request):
    user_groups = request.user.groups.values_list('name', flat=True)
    if user_groups and user_groups[0] == "samsung":
        return samsung_home(request)
    else:
        return generic_home(request, user_groups[0])


@login_required
def table(request):
    context = {'tags': get_sets_for_user(request.user, ["samsung"])}
    return render(request, 'table.html', context)


@login_required
def cards(request):
    context = {'tags': get_sets_for_user(request.user, ["samsung"])}
    return render(request, 'browse.html', context)


@login_required
def laser_table(request):
    context = {'tags': get_sets_for_user(request.user, ["laser"])}
    return render(request, 'table_laser.html', context)


@login_required
def laser_cards(request):
    context = {'tags': get_sets_for_user(request.user, ["laser"])}
    return render(request, 'browse_laser.html', context)


@login_required
def redox(request):
    # REMOVED CODE THAT Queries the DB for min and maxes.
    # context = RedoxPair.objects.all().aggregate(Min("redox_potential"),
    #                                             Max("redox_potential"),
    #                                             Max("substituent_count"),
    #                                             Max("water_solvation_energy"),
    #                                             Max("log_hyd_constant"))
    # for k in context.keys():
    #     if context[k] is not None:
    #         if k.endswith("max"):
    #             context[k] = round(float(context[k]), 2) + 0.01
    #         if k.endswith("min"):
    #             context[k] = round(float(context[k]), 2) - 0.01
    #     else:
    #         del(context[k])

    # set min max by hand
    # context = {"redox_potential__min": -1,
    #             "redox_potential__max": 3,
    #             "substituent_count__max": 8,
    #             "water_solvation_energy__max": 2,
    #             "water_solvation_energy__min": -5,
    #             "log_hyd_constant__max": 20,
    #             "log_hyd_constant__min": -20
    #             }

    context = {}
    context['tags'] = get_sets_for_user(request.user, ["flow_batt"])
    return render(request, 'browse_redox.html', context)


@login_required
def bubble(request):
    context = {'tags': get_sets_for_user(request.user, ["samsung"])}
    return render(request, 'browse.html', context)


@login_required
def redox_tree_json(request, redoxpair_id):
    initial_pair = get_object_or_404(RedoxPair, pk=redoxpair_id)
    family_tree = initial_pair.family_tree()
    return HttpResponse(json.dumps(family_tree), content_type="application/json")


@login_required
def redox_tree(request, redoxpair_id):
    context = {'redoxpair_id': redoxpair_id}
    return render(request, 'molvote/redox_tree.html', context)


@login_required
def redox_super_tree(request, count):
    all_pairs = RedoxPair.objects.filter(reduced__sets__name="0e")
    pairs = all_pairs.order_by('?')[:count]
    children = [pair.family_tree() for pair in pairs]
    context = {'children': json.dumps(children)}
    return render(request, 'molvote/redox_super_tree.html', context)


def details_helper(request, candidate):
    user_groups = request.user.groups.values_list('name', flat=True)
    if (not request.user.is_staff) and ((not candidate.released) or candidate.project not in user_groups):
        raise PermissionDenied
    comments = Comment.objects.filter(candidate=candidate)
    properties_dict = {}

    for prop in Property.objects.filter(candidate=candidate):
        method_list = properties_dict.get(prop.method.name, [])
        method_list.append(prop)
        properties_dict[prop.method.name] = method_list
    properties = [{'name': key, 'prop_list': value} for key, value in properties_dict.items()]

    context = {'candidate': candidate,
                'comments': comments,
                'properties': properties}
    return render(request, 'candidate.html', context)


@login_required
def detail(request, project_name, inchi_key):
    candidate = get_object_or_404(Candidate, project=project_name, inchi_key=inchi_key)
    return details_helper(request, candidate)


@login_required
def detail_by_id(request, candidate_id):
    candidate = get_object_or_404(Candidate, pk=candidate_id)
    return details_helper(request, candidate)


@login_required
def redox_graph_json(request, redoxpair_id):
    initial_pair = get_object_or_404(RedoxPair, pk=redoxpair_id)
    family_graph = initial_pair.family_graph()
    return HttpResponse(json.dumps(family_graph), content_type="application/json")


@login_required
def redox_graph(request, redoxpair_id):
    context = {'redoxpair_id': redoxpair_id}
    return render(request, 'molvote/redox_graph.html', context)


@login_required
def geometries(request, inchi_key):
    mol = mdm.Molecule.objects.filter(meta_data__inchi_key=inchi_key).first()
    conformers = [calc for calc in mol.calculation_list]
    all_calcs = []
    all_types = []
    if conformers:
        all_types.append(conformers[0].meta_data.worker_name)

    for conf in conformers:
        calcs_under_conf = {}
        calcs_to_check = [conf]
        while calcs_to_check:
            calc = calcs_to_check.pop(0)
            if calc.coord_list:
                calc_dict = {"name": calc.theory.theory_description,
                             "id": calc.id,
                             "geom": len(calc.coord_list) > 0}
                calcs_under_conf[calc.meta_data.worker_name] = calc_dict
                if calc.meta_data.worker_name not in all_types:
                    all_types.append(calc.meta_data.worker_name)
                for child in calc.child_calculation_list:
                    calcs_to_check.append(child)

        all_calcs.append(calcs_under_conf)

    if "pm7_opt" in all_types:
        all_types.remove("pm7_opt")
        all_types = [all_types[0], "pm7_opt"] + all_types[1:]
    output = []
    for ctype in all_types:
        calc_list = []
        for conf_dict in all_calcs:
            calc_list.append(conf_dict.get(ctype))
        output.append((ctype, calc_list))

    context = {'all_calcs': output, 'all_types': all_types}
    return render(request, 'molvote/geometries.html', context)


@only_staff_allowed
def statistics(request):
    announced = CandidateSet.objects.filter(announced=True, project="samsung")
    context = {'tags': [s.name for s in announced.order_by("name").only("name")]}
    return render(request, 'statistics.html', context)


def candidate_nickname_cmp(a, b):
    aname = a.candidate.nicknames.split(",")[-1]
    bname = b.candidate.nicknames.split(",")[-1]

    try:
        abatch = re.findall("[/a-z]+", aname)[0]
    except IndexError:
        abatch = ''
    try:
        bbatch = re.findall("[/a-z]+", bname)[0]
    except IndexError:
        bbatch = ''
    namecmp = cmp(abatch, bbatch)
    if namecmp != 0:
        return namecmp

    anums = [int(n) for n in re.findall("[0-9]+", aname)]
    bnums = [int(n) for n in re.findall("[0-9]+", bname)]
    try:
        namecmp = cmp(anums[0], bnums[0])
        if namecmp != 0:
            return namecmp
        return cmp(anums[1], bnums[1])
    except IndexError:
        return cmp(a.candidate.nicknames, b.candidate.nicknames)

nickname_key = functools.cmp_to_key(candidate_nickname_cmp)


@login_required
def ballot(request, ballot_id):
    ballot = get_object_or_404(Ballot, pk=ballot_id)

    context = {"tags": get_sets_for_user(request.user, ["samsung"]),
               "voting": True,
               "ballot_id": ballot.id,
               "ballot_name": ballot.name,
               "closed": ballot.closed(),
               "close_time": ballot.close_time,
               "tallies": count_votes(ballot)}
    return render(request, 'browse.html', context)


@login_required
def ballot_table(request, ballot_id):
    ballot = get_object_or_404(Ballot, pk=ballot_id)

    context = {"tags": get_sets_for_user(request.user, ["samsung"]),
               "voting": True,
               "ballot_id": ballot.id,
               "ballot_name": ballot.name,
               "closed": ballot.closed(),
               "close_time": ballot.close_time,
               "tallies": count_votes(ballot)}
    return render(request, 'table.html', context)


@login_required
def vote(request, ballot_id):
    rating, vote_id = request.POST["action"].split("-")
    vote = get_object_or_404(Vote, pk=vote_id)
    if vote.ballot.voter != request.user:
        raise PermissionDenied()
    ballot = get_object_or_404(Ballot, pk=ballot_id)
    if vote.ballot != ballot:
        raise PermissionDenied()
    if ballot.close_time is not None and timezone.now() > ballot.close_time:
        raise PermissionDenied("Ballot Closed")

    vote.rate(rating)

    response_data = {"candidate": vote.candidate.id,
                     "tallies": count_votes(ballot),
                     "id": vote.pk,
                     "rating": vote.rating}

    return HttpResponse(json.dumps(response_data), content_type="application/json")


class BallotManageView(StaffOnlyMixin, View):
    form_class = NewBallotForm
    template_name = 'ballot_admin.html'

    def ballot_dict(self, a):
        return {"name": a.name,
                "candidateset": a.candidateset,
                "voter": a.voter,
                "close_time": a.close_time,
                "count": Vote.objects.filter(ballot=a).count()}

    def get(self, request, **kwargs):
        form = self.form_class()
        all_ballots = Ballot.objects.all()
        ballots = [self.ballot_dict(a) for a in all_ballots]
        return render(request, self.template_name, {'form': form, 'all_ballots': ballots})

    def post(self, request, *args, **kwargs):
        form = self.form_class(request.POST)
        if form.is_valid():
            # <process form cleaned data>
            candidateset = form.cleaned_data['candidateset']
            name = form.cleaned_data['name']
            voter = form.cleaned_data['voter']
            close_time = form.cleaned_data['close_time']
            if Ballot.objects.filter(name=name, voter=voter).count() > 0:
                msg = "There is already a ballot with this name for {}"
                messages.error(request, msg.format(str(voter)))
            else:
                new_ballot = Ballot(name=name, voter=voter, close_time=close_time)
                new_ballot.save()
                new_ballot.import_candidateset(candidateset)
                messages.success(request, "New ballot successfully made for {}".format(str(voter)))

        all_ballots = Ballot.objects.all()
        ballots = [self.ballot_dict(a) for a in all_ballots]
        return render(request, self.template_name, {'form': form, 'all_ballots': ballots})


def cmp_up(a, b):
    up = cmp(a["votes"]["up"], b["votes"]["up"])
    if up != 0:
        return up
    else:
        return cmp(b["votes"]["down"], a["votes"]["down"])

up_key = functools.cmp_to_key(cmp_up)


class BallotResultsView(LoginRequiredMixin, View):
    template_name = 'ballot_results.html'

    def get(self, request, set_id):
        candidates = Candidate.objects.filter(sets__pk=set_id)
        cset = get_object_or_404(CandidateSet, pk=set_id)
        results = []
        ballot_info = {"voters": set(), "close_time": None}
        for c in candidates:
            votes = Vote.objects.filter(candidate=c)
            d = vars(c)
            vd = {"up": 0, "down": 0}
            for v in votes:
                if v.ballot:
                    voter_name = v.ballot.voter.username
                    ballot_info["voters"].add(voter_name)
                    if not ballot_info["close_time"]:
                        ballot_info["close_time"] = v.ballot.close_time
                    vd[voter_name] = v.rating
                    if v.rating is not None and v.rating > 0:
                        vd["up"] += 1
                    if v.rating is not None and v.rating < 0:
                        vd["down"] += 1
            d["short_inchi"] = d["inchi_key"].split("-")[0]
            d["votes"] = vd
            results.append(d)

        context = {'ballot_info': ballot_info,
                   'cset': cset,
                   'candidates': reversed(sorted(results, key=up_key))}
        return render(request, self.template_name, context)
