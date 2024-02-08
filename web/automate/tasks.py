from __future__ import absolute_import
import datetime
import os
import time
import logging
from logging.handlers import TimedRotatingFileHandler

from django.core.management import call_command
from django.core.cache import cache
import django.utils.timezone
from django.conf import settings
from django.db import DatabaseError
from django.core.exceptions import ObjectDoesNotExist
from celery import task
from celery.utils.log import get_task_logger
from mongoengine.context_managers import switch_db

from tools import extract_results_from_db as extract
from tools import extract_quinone_results_from_db as extract_quinones
from tools import load_quinones
from tools import update_quinones_sql
from tools.propagate import Propagation

from aag_python.molecular_storage import molecular_data_models as mdm
from molvote.models import CandidateSet, Candidate, Method, Property


logger = get_task_logger(__name__)

LOCK_EXPIRE = 3600 * 24  # Lock expires in 24 hours

## Use django database config to associate different projects with mongo databases
DB_ALIAS_LOOKUP = {}
for db_config in settings.DATABASES.values():
    for project_name in db_config.get("PROJECTS", []):
        DB_ALIAS_LOOKUP[project_name] = db_config["ALIAS"]

def setup_log(name, path, streaming):
    logger = logging.getLogger(name)
    if not len(logger.handlers):
        logger.setLevel(logging.WARNING)
        # create file handler which logs even debug messages

        filehandler = TimedRotatingFileHandler(os.path.expanduser(path), 'midnight', 1, backupCount=3)
        filehandler.suffix = "-%y-%m-%d"
        logger.addHandler(filehandler)

        filehandler.setLevel(logging.ERROR)
        # create console handler with a higher log level
        streamhandler = logging.StreamHandler()
        streamhandler.setLevel(logging.ERROR)
        # create formatter and add it to the handlers
        formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M')
        streamhandler.setFormatter(formatter)
        filehandler.setFormatter(formatter)
        # add the handlers to logger
        if streaming:
            logger.addHandler(streamhandler)
        logger.addHandler(filehandler)
    return logger

propagate_logger = setup_log('propagate', '/var/log/celery/propagate.log', False)


@task(name="automate.tasks.dump_csv")
def dump_csv():
    extract.dump_excited_state_w_lifetime(mdm.Theory.TheoryLevelValue.TDDFT_HYBRID_B3LYP,
                                          path=settings.DUMP_CSV_PATH)

USE_KEYS = ("inchi_key",
            "smiles",
            "nicknames",
            "date",
            "sascore",
            "weight",
            )

PROPERTIES = ["absorption",
              "splitting",
              "strength",
              "rate",
              "homo",
              "lumo",
              "total_energy",
              "rmsd_s0_t1",
              "reorg_en_t1",
              "reorg_en_s0",
              "s1_rate"]

CALC_TRANSFER_MAP_LIST = [  {"method": "b3lyp_tddft_631gs_rpa_s0_geom",
                             "parents": ['b3lyp_6-31gs_opt', 'b3lyp_6-31gs_opt_qchem'],
                             "worker_name": "b3lyp_6-31gs_tddft",
                             "level": "tddft_hyb_b3l"},
                            {"method": "b3lyp_tddft_631gs_rpa_t1_geom",
                             "parents": ['b3lyp_6-31gs_opt_t1_qchem', 'b3lyp_6-31gs_opt_t1'],
                             "worker_name": "b3lyp_6-31gs_tddft",
                             "level": "tddft_hyb_b3l"},
                            {"method": "b3lyp_tddft_631gs_pcm_rpa_s0_geom",
                             "parents": ['b3lyp_6-31gs_opt', 'b3lyp_6-31gs_opt_qchem'],
                             "worker_name": "b3lyp_6-31gs_tddft_pcm",
                             "level": "tddft_hyb_b3l"},
                            {"method": "b3lyp_tddft_631gs_pcm_rpa_t1_geom",
                             "parents": ['b3lyp_6-31gs_opt_t1_qchem', 'b3lyp_6-31gs_opt_t1'],
                             "worker_name": "b3lyp_6-31gs_tddft_pcm",
                             "level": "tddft_hyb_b3l"},
                            {"method": "camb3lyp_tddft_631gs_rpa_s0_geom",
                             "parents": ['b3lyp_6-31gs_opt', 'b3lyp_6-31gs_opt_qchem'],
                             "worker_name": "camb3lyp_6-31gs_tddft",
                             "level": "tddft_hyb_camb3l"},
                            {"method": "camb3lyp_tddft_631gs_rpa_t1_geom",
                             "parents": ['b3lyp_6-31gs_opt_t1_qchem', 'b3lyp_6-31gs_opt_t1'],
                             "worker_name": "camb3lyp_6-31gs_tddft",
                             "level": "tddft_hyb_camb3l"},
                            {"method": "wb97xd_tddft_631gs_rpa_s0_geom",
                             "parents": ['b3lyp_6-31gs_opt', 'b3lyp_6-31gs_opt_qchem'],
                             "worker_name": "wb97xd_6-31gs_tddft",
                             "level": "tddft_hyb_wbx"},
                            {"method": "wb97xd_tddft_631gs_rpa_t1_geom",
                             "parents": ['b3lyp_6-31gs_opt_t1_qchem', 'b3lyp_6-31gs_opt_t1'],
                             "worker_name": "wb97xd_6-31gs_tddft",
                             "level": "tddft_hyb_wbx"},
                            {"method": "m062x_tddft_631gs_rpa_s0_geom",
                             "parents": ['b3lyp_6-31gs_opt', 'b3lyp_6-31gs_opt_qchem'],
                             "worker_name": "m062x_6-31gs_tddft",
                             "level": "tddft_hyb_m062x"},
                            {"method": "m062x_tddft_631gs_rpa_t1_geom",
                             "parents": ['b3lyp_6-31gs_opt_t1_qchem', 'b3lyp_6-31gs_opt_t1'],
                             "worker_name": "m062x_6-31gs_tddft",
                             "level": "tddft_hyb_m062x"},
                            {"method": "lcwpbeh_tddft_631gs_rpa_s0_geom",
                             "parents": ['b3lyp_6-31gs_opt', 'b3lyp_6-31gs_opt_qchem'],
                             "worker_name": "lcwpbeh_6-31gs_tddft",
                             "level": "tddft_lrc_hyb_lcwpbeh"},
                            {"method": "lcwpbeh_tddft_631gs_rpa_t1_geom",
                             "parents": ['b3lyp_6-31gs_opt_t1_qchem', 'b3lyp_6-31gs_opt_t1'],
                             "worker_name": "lcwpbeh_6-31gs_tddft",
                             "level": "tddft_lrc_hyb_lcwpbeh"},
                            {"method": "bhhlyp_tddft_631gs_rpa_s0_geom",
                             "parents": ['b3lyp_6-31gs_opt', 'b3lyp_6-31gs_opt_qchem'],
                             "worker_name": "bhhlyp_6-31gs_tddft",
                             "level": "tddft_hyb_bhh"},
                            {"method": "bhhlyp_tddft_631gs_rpa_t1_geom",
                             "parents": ['b3lyp_6-31gs_opt_t1_qchem', 'b3lyp_6-31gs_opt_t1'],
                             "worker_name": "bhhlyp_6-31gs_tddft",
                             "level": "tddft_hyb_bhh"},
                            {"method": "b3lyp_tddft_631gs_rpa_s1_geom",
                             "parents": ['b3lyp_6-31gs_opt_s1'],
                             "worker_name": "b3lyp_6-31gs_tddft",
                             "level": "tddft_hyb_b3l"},
                            {"method": "b3lyp_tddft_631gs_pcm_rpa_s1_geom",
                             "parents": ['b3lyp_6-31gs_opt_s1'],
                             "worker_name": "b3lyp_6-31gs_tddft_pcm",
                             "level": "tddft_hyb_b3l"},
                            {"method": "camb3lyp_tddft_631gs_rpa_s1_geom",
                             "parents": ['b3lyp_6-31gs_opt_s1'],
                             "worker_name": "camb3lyp_6-31gs_tddft",
                             "level": "tddft_hyb_camb3l"},
                            {"method": "wb97xd_tddft_631gs_rpa_s1_geom",
                             "parents": ['b3lyp_6-31gs_opt_s1'],
                             "worker_name": "wb97xd_6-31gs_tddft",
                             "level": "tddft_hyb_wbx"},
                            {"method": "m062x_tddft_631gs_rpa_s1_geom",
                             "parents": ['b3lyp_6-31gs_opt_s1'],
                             "worker_name": "m062x_6-31gs_tddft",
                             "level": "tddft_hyb_m062x"},
                            {"method": "lcwpbeh_tddft_631gs_rpa_s1_geom",
                             "parents": ['b3lyp_6-31gs_opt_s1'],
                             "worker_name": "lcwpbeh_6-31gs_tddft",
                             "level": "tddft_lrc_hyb_lcwpbeh"},
                            {"method": "bhhlyp_tddft_631gs_rpa_s1_geom",
                             "parents": ['b3lyp_6-31gs_opt_s1'],
                             "worker_name": "bhhlyp_6-31gs_tddft",
                             "level": "tddft_hyb_bhh"}
                            ]


def create_or_update_candidate(data, project, method_name):
    args = dict((k, data[k]) for k in USE_KEYS)
    # patch for change in Candidate fields
    if "date" in args:
        aware_time = django.utils.timezone.make_aware(django.utils.dateparse.parse_datetime(args["date"]),
                                                      django.utils.timezone.utc)
        args["calc_time"] = aware_time
        del(args["date"])
    args["project"] = project

    # args_for_new_candidate = args.copy()
    # args_for_new_candidate.update(dict((k, data[k]) for k in PROPERTIES))
    method, created = Method.objects.get_or_create(name=method_name)
    cquery = Candidate.objects.filter(inchi_key=data["inchi_key"], project=project)
    if len(cquery) > 1:
        msg = "Found more than one Candidate for key {} and project {}".format(data["inchi_key"], project)
        raise StandardError(msg)
    if not cquery.exists():
        can = Candidate(**args)
        can.save()
        candidate_created = True
    else:
        cquery.update(**args)
        can = cquery[0]
        candidate_created = False

    for prop_name in PROPERTIES:
        if prop_name in data and data[prop_name] is not None:
            try:
                prop = Property.objects.get(method=method, name=prop_name, candidate=can)
            except ObjectDoesNotExist:
                prop = Property(method=method, name=prop_name, candidate=can)
            prop.value = data[prop_name]
            prop.save()

    # handle sets
    set_names = [s.name for s in can.sets.all()]
    tag_names = data["tags"].split(" ")
    for tagname in tag_names:
        if tagname not in set_names:
            cset, set_created = CandidateSet.objects.get_or_create(name=tagname, project=project)
            can.sets.add(cset)

    for cset in can.sets.all():
        if cset.name not in tag_names:
            can.sets.remove(cset)

    donor_and_acceptor_found = True
    try:
        can.donor = Candidate.objects.get(inchi_key=data["donor_inchi_key"], project=project)
    except (KeyError, Candidate.DoesNotExist):
        donor_and_acceptor_found = False
    try:
        can.acceptor = Candidate.objects.get(inchi_key=data["acceptor_inchi_key"], project=project)
    except (KeyError, Candidate.DoesNotExist):
        donor_and_acceptor_found = False
    try:
        can.bridge1 = Candidate.objects.get(inchi_key=data["bridge1_inchi_key"], project=project)
    except (KeyError, Candidate.DoesNotExist):
        pass
    try:
        can.bridge2 = Candidate.objects.get(inchi_key=data["bridge2_inchi_key"], project=project)
    except (KeyError, Candidate.DoesNotExist):
        pass

    saved = False
    try_count = 0
    while not saved:
        try:
            can.save()
            saved = True
        except DatabaseError:
            try_count += 1
            if try_count > 10:
                raise
            time.sleep(0.01)

    return candidate_created, donor_and_acceptor_found


@task(name="automate.tasks.update_sql")
def update_sql(inchis=[],
               tag_list=[],
               src_project="samsung",
               target_project="samsung",
               within_last_hours=None,
               include_transferred=False,
               mark_as_transferred=True,
               ignore_lock=False):
    lock_id = 'update_sql_lock-{}'.format(target_project)

    # cache.add fails if if the key already exists

    def acquire_lock():
        return cache.add(lock_id, 'true', LOCK_EXPIRE)
    # delete is very slow, but we have to use it to take
    # advantage of using add() for atomic locking

    def release_lock():
        cache.delete(lock_id)

    if not ignore_lock and not acquire_lock():
        logger.warning("last task with lock [{}] unfinished.  Skipping work round.".format(lock_id))
        return
    try:
        if within_last_hours:
            from_datetime = datetime.datetime.utcnow() - datetime.timedelta(hours=within_last_hours)
        else:
            from_datetime = None

        logger.info("Will move mols since {}".format(from_datetime))
        created = 0
        updated = 0

        for map_dict in CALC_TRANSFER_MAP_LIST:
            query = extract.query_tddft_calcs(inchis=inchis,
                                              level=map_dict["level"],
                                              tag_list=tag_list,
                                              from_datetime=from_datetime,
                                              with_ancestry=True,
                                              project=src_project,
                                              include_transferred=include_transferred,
                                              worker_name=map_dict["worker_name"],
                                              logger=logger,
                                              parent_list=map_dict["parents"])
            for data in query:
                cand_created, parents_found = create_or_update_candidate(data, target_project, map_dict["method"])
                if cand_created:
                    created += 1
                else:
                    updated += 1
                if mark_as_transferred:
                    mdm.Calculation.objects.filter(id__in=data["ids"]).update(set__meta_data__status__transferred=True)
        msg = "Candidate objects. New: {}. Updated: {}".format(created, updated)
        logger.info(msg)
        print msg

    finally:
        if not ignore_lock:
            release_lock()


@task(name="automate.tasks.update_flow_batt")
def update_flow_batt(inchis=[], tag_list=[],
                     db_alias="flow_batt",
                     override_smiles_confirm=False,
                     minimums_only=True,
                     family='quinone',
                     ignore_lock=False,
                     michael_mode='only_H'):
    lock_id = 'update_sql_lock-flow_batt'
    # cache.add fails if if the key already exists

    def acquire_lock():
        return cache.add(lock_id, 'true', LOCK_EXPIRE)
    # delete is very slow, but we have to use it to take
    # advantage of using add() for atomic locking

    def release_lock():
        cache.delete(lock_id)

    if not ignore_lock and not acquire_lock():
        logger.warning("last task with lock [{}] unfinished.  Skipping work round.".format(lock_id))
        return
    try:
        with switch_db(mdm.Calculation, db_alias):
            with switch_db(mdm.Molecule, db_alias):
                for redox_dict in extract_quinones.query_redox_calcs(inchis, tag_list,
                                                                     minimums_only=minimums_only,
                                                                     override_smiles_confirm=override_smiles_confirm,
                                                                     family=family,
                                                                     michael_mode='only_H'):
                    load_quinones.add_redox(redox_dict)
        tidier = update_quinones_sql.QuinoneSqlTidyUp(inchis=inchis, sets=tag_list)
        tidier.query_quinones()
        tidier.set_minima()
        tidier.tag_n_subs()
        tidier.make_pair_of_pairs()
    finally:
        release_lock()


@task(name="backup_databases")
def backup_databases():
    call_command("dbbackup")


@task(name="propagate_oled")
def propagate_oled(projects=['samsung'], ignore_lock=False):
    lock_id = 'propagate_sql_lock-{}'.format(projects)

    # cache.add fails if if the key already exists

    def acquire_lock():
        return cache.add(lock_id, 'true', LOCK_EXPIRE)

    # delete is very slow, but we have to use it to take
    # advantage of using add() for atomic locking

    def release_lock():
        cache.delete(lock_id)

    if not ignore_lock and not acquire_lock():
        logger.warning("last task with lock [{}] unfinished.  Skipping work round.".format(lock_id))
        return

    try:
        Propagation(logger=propagate_logger).run_propagator()
    finally:
        release_lock()


@task(name="propagate_flow_batt")
def propagate_flow_batt(projects=['flow_batt'], ignore_lock=False):
    lock_id = 'propagate_sql_lock-{}'.format(projects)

    # cache.add fails if if the key already exists

    def acquire_lock():
        return cache.add(lock_id, 'true', LOCK_EXPIRE)

    # delete is very slow, but we have to use it to take
    # advantage of using add() for atomic locking

    def release_lock():
        cache.delete(lock_id)

    if not ignore_lock and not acquire_lock():
        logger.warning("last task with lock [{}] unfinished.  Skipping work round.".format(lock_id))
        return

    try:
        for proj in projects:
            with switch_db(mdm.Calculation, DB_ALIAS_LOOKUP[proj]):
                with switch_db(mdm.Molecule, DB_ALIAS_LOOKUP[proj]):
                    Propagation(logger=propagate_logger).run_propagator(project=proj)
    finally:
        release_lock()
