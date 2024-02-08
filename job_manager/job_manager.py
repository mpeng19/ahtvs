from __future__ import print_function

"""
Job Manager

Tim Hirzel <hirzel@chemistry.harvard.edu>
April 2014
A poorly named tool that watches a inbox directory,
automatically adds job scripts to a queue,
(ie. slurm) monitors their progress, and finally
moves job results into a completed directory.

More info in README.
"""
import argparse
import datetime
import importlib
import logging
import glob
import json
import os
import re
import shutil
import signal
import stat
import tarfile
import time

from signalhandler import SignalHandler

INBOX_NAME = "inbox"
PENDING_NAME = "pending"
COMPLETED_NAME = "completed"

DEFAULT_MAX_QUEUE_LENGTH = 100
DEFAULT_MAX_PENDING_LENGTH = 20
DEFAULT_PERIOD_SECONDS = 120
DEFAULT_MAX_ERROR_COUNT = 0
DEFAULT_SUBMIT_DELAY = 1
SUBMISSION_PROBLEM_DELAY_MULTIPLE = 10  # when submission fails. mutiply submission rate by this
DEFAULT_SUBMIT_DELAY_LOCAL = 0.05
DEFAULT_JOB_ID_FILENAME = "job_manager-job_id"
DEFAULT_JOB_FILENAME = "job.sh"
DEFAULT_LOG_FILENAME = "job_manager-log"
DEFAULT_ERROR_FILENAME = "job_manager-errors"
DEFAULT_JOB_NAME_FILENAME = "job_name.txt"
DEFAULT_DO_CPU_GPU = False
DEFAULT_COMPLETE_FILENAME = "job_manager-complete"
DEFAULT_FINISHED_CHECKPOINT_FILENAME = "job_manager-finished"
DEFAULT_ENGINE_MIN_SCAN_INTERVAL = 5 * 60 # 5 minutes

DEFAULT_LOG_LEVEL = logging.INFO

JOB_CONFIRMED_STR = "confirmed"

APPEND_LINES_TO_JOB = """

# LINES BELOW ARE ADDED AUTOMATICALLY BY JOB MANAGER #
touch {}

"""


class JobManager:

    def __init__(
        self,
        root_path,
        engine,
        period=DEFAULT_PERIOD_SECONDS,
        submit_delay=DEFAULT_SUBMIT_DELAY,
        max_queue_len=DEFAULT_MAX_QUEUE_LENGTH,
        max_pending_len=DEFAULT_MAX_PENDING_LENGTH,
        max_error_count=DEFAULT_MAX_ERROR_COUNT,
        setup=False,
        job_filename=DEFAULT_JOB_FILENAME,
        complete_filename=DEFAULT_COMPLETE_FILENAME,
        error_filename=DEFAULT_ERROR_FILENAME,
        job_id_filename=DEFAULT_JOB_ID_FILENAME,
        job_name_filename=DEFAULT_JOB_NAME_FILENAME,
        do_cpu_gpu=DEFAULT_DO_CPU_GPU,
        output_group_writable=False,
        presubmit_checkpoint_filename=None,
        finished_checkpoint_filename=None,
        submission_grace_period=None,
        global_lock_file=None,
        max_error_percent=None,
        error_lookback_time=None,
        error_cooldown_time=None,
        error_lookback_thresh=None,
        scratch_dir=None,
    ):
        self.root_path = root_path
        self.engine = engine
        self.job_filename = job_filename
        self.period = period
        self.do_cpu_gpu = do_cpu_gpu
        self.default_submit_delay = submit_delay
        self.current_submit_delay = submit_delay
        self.max_queue_len = max_queue_len
        self.max_pending_len = max_pending_len
        self.max_error_count = max_error_count
        self.job_id_filename = job_id_filename
        self.presubmit_checkpoint_filename = presubmit_checkpoint_filename
        self.finished_checkpoint_filename = finished_checkpoint_filename
        self.complete_filename = complete_filename
        self.error_filename = error_filename
        self.job_name_filename = job_name_filename
        self.output_group_writable = output_group_writable
        self.submission_grace_period = submission_grace_period
        self.global_lock_file = global_lock_file
        self.max_error_percent = max_error_percent
        self.error_lookback_time = error_lookback_time
        self.error_cooldown_time = error_cooldown_time
        self.error_lookback_thresh = error_lookback_thresh

        self.scratch_dir = scratch_dir
        os.makedirs(self.scratch_dir, exist_ok=True)

        self.queued_jobs = {}
        self.recent_submissions = {}
        self.recent_completions = {}
        self.inbox_size = 0
        self.cooling_down = False
        self.tick_counter = 0
        self.tar_info_name = 'tar_info.json'

        dirs = os.listdir(self.root_path)
        expect = [INBOX_NAME, PENDING_NAME, COMPLETED_NAME]
        for dname in expect:
            if dname not in dirs:
                if setup:
                    os.mkdir(os.path.join(self.root_path, dname))
                else:
                    err = "directory {0} not found".format(os.path.join(self.root_path, dname))
                    self._log(err, level=logging.CRITICAL)
                    raise Exception(err)

        self.cpu_gpu_check_period = 0
        self.inbox_dir = os.path.abspath(os.path.join(self.root_path, INBOX_NAME))
        self.pending_dir = os.path.abspath(os.path.join(self.root_path, PENDING_NAME))
        self.completed_dir = os.path.abspath(os.path.join(self.root_path, COMPLETED_NAME))

    def submit_from_inbox(self):
        """
        submit job dirs from inbox to pending dir and start jobs
        """

        tar_completions = self.process_tars_in_inbox()
        self.recent_completions.update(tar_completions)

        just_submitted_ids = {}
        job_dirs = [d for d in os.listdir(self.inbox_dir) if os.path.isdir(os.path.join(self.inbox_dir, d))]
        job_dirs.sort()
        self.inbox_size = len(job_dirs)

        total_len, pending_len = self.engine.get_pending_count(pending_only=True)
        open_spots = min(self.max_queue_len - total_len, self.max_pending_len - pending_len)
        open_spots = max(open_spots, 0)
        self._log("{0} open spots in queue".format(open_spots))
        for jobdir in job_dirs[:open_spots]:
            if self.presubmit_checkpoint_filename:
                presubmit_checkpoint_path = os.path.join(
                    self.inbox_dir, jobdir, self.presubmit_checkpoint_filename)
                if not os.path.isfile(presubmit_checkpoint_path):
                    self._log(
                        "presubmit checkpoint file '{}' not found".format(
                            presubmit_checkpoint_path), level=logging.DEBUG)
                    continue
            job_script_path = os.path.join(self.inbox_dir, jobdir, self.job_filename)
            if not os.path.isfile(job_script_path):
                self._log(
                    "Could not find job script: {0}".format(job_script_path),
                    level=logging.DEBUG)
                continue
            if not os.access(os.path.join(self.inbox_dir, jobdir), os.W_OK):
                self._log(
                    "No write access to directory: {}".format(
                        os.path.join(self.inbox_dir, jobdir)),
                    level=logging.WARNING)
                continue
            if not os.access(job_script_path, os.W_OK):
                self._log("No write access to file: {}".format(
                    job_script_path), level=logging.WARNING)
                continue
            pending_path = os.path.join(self.pending_dir, jobdir)
            if os.path.exists(pending_path):
                self._log(
                    "Already job dir of same name in pending: {}".format(
                        pending_path), level=logging.ERROR)
                continue
            with SignalHandler([signal.SIGINT]) as handler:
                shutil.move(os.path.join(self.inbox_dir, jobdir), pending_path)
                job_script_path = os.path.join(pending_path, self.job_filename)
                with open(job_script_path, mode='a') as f:
                    f.write(APPEND_LINES_TO_JOB.format(self.complete_filename))
                job_name_file_path = os.path.join(pending_path, self.job_name_filename)
                if os.path.isfile(job_name_file_path):
                    try:
                        with open(job_name_file_path, mode='r') as name_file:
                            jobname = ["_"].join(name_file.readlines()[0].split())
                    except Exception as e:
                        self._log(
                            "Unable to read job name file: {0}".format(e),
                            level=logging.WARNING)
                        jobname = jobdir
                else:
                    jobname = jobdir
                jid = self.engine.submit_job(job_script_path, name=jobname)
                if jid is None or not re.match('\d+', jid.strip()):
                    error_file_path = os.path.join(pending_path, self.error_filename)
                    with open(error_file_path, 'w') as f:
                        f.write("Job submission failed, no job id returned")
                    self._move_from_pending_to_completed(jobdir)
                    self.current_submit_delay = self.default_submit_delay * SUBMISSION_PROBLEM_DELAY_MULTIPLE
                    msg = "Failed starting job in {0}. Slowing submission delay to {1} seconds"
                    msg = msg.format(jobdir, self.current_submit_delay)
                    self._log(msg, level=logging.ERROR)
                else:
                    jid = jid.strip()
                    jid_file_path = os.path.join(pending_path, self.job_id_filename)
                    with open(jid_file_path, 'w') as f:
                        f.write(jid)
                    just_submitted_ids[jid] = time.time()
                    self._log("Job Submitted: {0}".format(jobdir))
                    self.current_submit_delay = self.default_submit_delay
                if handler.interrupted:
                    self._log(
                        "Terminated with signal {0}".format(handler.sig_used),
                        level=logging.FATAL)
                    break

            time.sleep(self.current_submit_delay)
        return just_submitted_ids

    def process_tars_in_inbox(self):
        completions = {}
        tars = self.list_tars(dir=self.inbox_dir)
        for tar in tars:
            tar_meta = self.get_tar_meta(tar=tar)
            try:
                self.validate_job_tar(tar=tar, tar_meta=tar_meta)
                self.unpack_tar_to_inbox_job_dir(tar=tar, tar_meta=tar_meta)
            except Exception as error:
                self._log(
                    """Could not process tarfile '{tar}', moving to
                    completed. Error was {error}""".format(tar=tar, error=error)
                )
                dest_tar = os.path.join(
                    self.completed_dir, tar_meta['tar_name'])
                shutil.move(tar, dest_tar)
                completions[tar] = {
                    'status': 'FAILED',
                    'time': time.time()
                }
        return completions

    def list_tars(self, dir=None):
        return glob.glob(dir + '/*.tar.*z')

    def get_tar_meta(self, tar=None):
        tar_name = os.path.basename(tar)
        base_name = re.match(r'(.*)\.tar(\..*)?', tar_name).group(1)
        compression_type = tar_name.split('.')[-1]
        tar_meta = {
            'original_submission_path': tar,
            'tar_name': tar_name,
            'base_name': base_name,
            'compression_type': compression_type,
        }
        return tar_meta

    def validate_job_tar(self, tar=None, tar_meta=None):
        try:
            handle = tarfile.open(name=tar,
                                  mode='r:%s' % tar_meta['compression_type'])
            handle.getmember(self.job_filename)
        except KeyError:
            raise Exception("Unable to find file '{job_filename}'".format(
                self.job_filename))

    def unpack_tar_to_inbox_job_dir(self, tar=None, tar_meta=None):
        dest_dir = os.path.join(self.inbox_dir, tar_meta['base_name'])
        if os.path.exists(dest_dir):
            msg = ("""Ignoring tarfile {tar}, because it would unpack to
                   existing dir '{dest_dir}'""").format(
                       tar=tar,
                       dest_dir=dest_dir
                   )
            self._log(msg)
            return
        timestamp = time.time()
        tmp_tar = os.path.join(self.scratch_dir, tar_meta['tar_name'])
        shutil.move(tar_meta['original_submission_path'], tmp_tar) 
        tmp_dir = "{tmp_dir_path}.{timestamp}".format(
            tmp_dir_path=os.path.join(self.scratch_dir, tar_meta['base_name']),
            timestamp=timestamp,
        )
        os.makedirs(tmp_dir)
        handle = tarfile.open(name=tmp_tar,
                              mode='r:%s' % tar_meta['compression_type'])
        handle.extractall(path=tmp_dir)
        tar_info_path = os.path.join(tmp_dir, self.tar_info_name)
        self.write_tar_info(output_path=tar_info_path, tar_meta=tar_meta)
        if self.presubmit_checkpoint_filename:
            checkpoint_path = os.path.join(
                tmp_dir, self.presubmit_checkpoint_filename)
            open(checkpoint_path, 'a').close()
        shutil.move(tmp_dir, dest_dir)

    def write_tar_info(self, output_path=None, tar_meta=None):
        tar_info = {k: v for k, v in tar_meta.items() if k != 'handle'}
        with open(output_path, 'w') as f: json.dump(tar_info, f)

    def _pending_job_id(self, jobdir):
        job_id_path = os.path.join(self.pending_dir, jobdir, self.job_id_filename)
        with open(job_id_path, 'r') as fp:
            lines = fp.readlines()
        try:
            jid = lines[0].strip()
            if jid == "":
                raise Exception("Job found with no ID")
            else:
                return jid
        except IndexError:
            raise Exception("Job found with no ID")

    def check_for_completions(self):
        """
        move any finished jobs from pending dir to completed dir.
        """
        completions = {}
        for jobdir in os.listdir(self.pending_dir):
            full_dir = os.path.join(self.pending_dir, jobdir)
            job_id_path = os.path.join(full_dir, self.job_id_filename)
            try:
                with open(job_id_path, 'r') as f: jid = f.readlines()[0].strip()
            except Exception as e:
                self._log("Could not read job_id file '{}'".format(
                   job_id_path))
                continue
            if (jid in self.recent_submissions) or (jid in self.queued_jobs):
                continue
            complete_checkpoint_path = os.path.join(full_dir,
                                                    self.complete_filename)
            if not os.path.isfile(complete_checkpoint_path):
                self._log(
                    """Could not find completion checkpoint file '{}', job 
                    may have failed. Moving jobdir to completed""".format(
                        complete_checkpoint_path))
                self._move_from_pending_to_completed(jobdir)
                completion_status = 'FAILED'
            else:
                self._log("Job Completed: {}".format(jobdir))
                self._move_from_pending_to_completed(jobdir)
                completion_status = 'SUCCEEDED'
            completions[jid] = {
                'status': completion_status,
                'time': time.time()
            }
        return completions

    def _log(self, msg, level=logging.INFO):
        logging.log(level, msg)
        print("\n" + msg)

    def _move_from_pending_to_completed(self, jobdir):
        source_dir = os.path.join(self.pending_dir, jobdir)
        tar_info_path = os.path.join(source_dir, self.tar_info_name)
        if os.path.exists(tar_info_path):
            tar_info = json.load(open(tar_info_path))
            dest_tar = os.path.join(self.completed_dir, tar_info['tar_name'])
            self._move_dir_to_tar(source_dir=source_dir, tar_info=tar_info,
                                  dest_tar=dest_tar)
        else:
            dest_dir = os.path.join(self.completed_dir, jobdir)
            self._move_dir(
                source_dir=source_dir,
                dest_dir=dest_dir,
                checkpoint_filename=self.finished_checkpoint_filename
            )

    def _move_dir_to_tar(self, source_dir=None, tar_info=None, dest_tar=None):
        tmp_tar = os.path.join(self.scratch_dir, tar_info['tar_name'])
        tar_handle = tarfile.open(
            name=tmp_tar,
            mode='w:{}'.format(tar_info['compression_type'])
        )
        tar_handle.add(source_dir, arcname='')
        tar_handle.close()
        shutil.move(tmp_tar, dest_tar)
        shutil.rmtree(source_dir)

    def _move_dir(self, source_dir=None, dest_dir=None,
                  checkpoint_filename=None):
        with SignalHandler([signal.SIGINT]) as handler:
            if self.output_group_writable:
                stat_info = os.stat(source_dir)
                gid = stat_info.st_gid
                my_uid = os.geteuid()
                for dir_root, dirs, files in os.walk(source_dir):
                    for f in dirs + files:
                        path = os.path.join(dir_root, f)
                        stat_info = os.stat(path)
                        if stat_info.st_uid == my_uid:
                            os.chmod(path, stat_info.st_mode | stat.S_IWGRP | stat.S_IWOTH)
                            os.chown(path, stat_info.st_uid, gid)
            shutil.move(source_dir, dest_dir)

            if handler.interrupted:
                msg = "Terminated with signal {0}".format(handler.sig_used)
                logging.warning(msg)
                raise Exception(msg)
            else:
                if checkpoint_filename:
                    checkpoint_path = os.path.join(dest_dir,
                                                   checkpoint_filename)
                    open(checkpoint_path, 'a').close()

    def check_cpu_in_gpu(self, cpu_part, gpu_part, all_parts, qos, max_jobs, nsec, exc_nodes=None):
        if engine_class_name != "SlurmEngine":
            raise AssertionError('cpu_to_gpu only allowed at Odyssey')

        gpu_run_list = self.engine._check_queue('RUNNING', partition=gpu_part, qos=None, all_users=True)
        self._log('There are {0} jobs running in the gpu \
                            queues'.format(len(gpu_run_list)))
        cpu_list = self.engine._check_queue('PENDING', cpu_part, qos='normal')
        self._log('There are {0} jobs pending in the cpu \
                            queues'.format(len(cpu_list)))
        cpu_list.sort(key=lambda x: x[2])
        if len(gpu_run_list) < max_jobs and (time.time() - self.cpu_gpu_check_period) > nsec:
            self.cpu_gpu_check_period = time.time()
            cpu_jobs_to_swap = min(len(cpu_list), (max_jobs-len(gpu_run_list))/2)
            self._log("Found {0} free cpu spots in gpu partitions; \
                                moving {1} jobs  to {2}".format(
                                    (max_jobs - len(gpu_run_list)),
                                    cpu_jobs_to_swap, gpu_part))
            for i in range(cpu_jobs_to_swap):
                self.engine.update_job(cpu_list[i][1], all_parts, qos, exc_node_list_string=exc_nodes)
                time.sleep(self.current_submit_delay)
                self._log("Job {0} swapped to partition {1}".format(
                    cpu_list[i][1], gpu_part))

    def has_global_lock(self):
        has_lock = (self.global_lock_file and \
                    os.path.isfile(self.global_lock_file))
        return has_lock

    def prune_recent_submissions(self):
        now = time.time()
        keys_to_remove = set()
        for jid, registration_time in list(self.recent_submissions.items()):
            exceeds_submission_grace_period = (now - registration_time) > \
                    self.submission_grace_period
            if exceeds_submission_grace_period: keys_to_remove.add(jid)
        self.recent_submissions = {
            k: v for k, v in self.recent_submissions.items()
            if k not in keys_to_remove
        }

    def prune_recent_completions(self):
        now = time.time()
        keys_to_remove = set()
        for job_id, completion in self.recent_completions.items():
            outside_lookback_window = (now - completion['time'] ) > \
                    self.error_lookback_time
            if outside_lookback_window: keys_to_remove.add(job_id)
        self.recent_completions = {
            k: v for k, v in self.recent_completions.items()
            if k not in keys_to_remove
        }

    def get_error_percent(self):
        if len(self.recent_completions) == 0:
            error_percent = 0
        else:
            num_recent_errors = self.get_num_recent_errors()
            error_percent = 100.0 * \
                    (num_recent_errors / len(self.recent_completions))
        return error_percent

    def get_num_recent_errors(self):
        num_errors = 0
        for completion in self.recent_completions.values():
            if completion['status'] == 'FAILED': num_errors += 1
        return num_errors

    def update_queued_jobs(self):
        queued_job_list = self.engine.get_pending()
        self.queued_jobs = dict((str(j[1]), j) for j in queued_job_list)

    def run(self):
        self.tick_counter = 0
        while True:
            self.tick_counter += 1
            self._log(
                "{oomph} TICK #{tick_counter}, {time} {oomph}".format(
                    oomph=('=' * 10),
                    tick_counter=self.tick_counter,
                    time=datetime.datetime.now()))
            if self.has_global_lock():
                self._log(
                    "skipping tick, lock file '{lock_file}' exists.".format(
                        lock_file=self.global_lock_file))
            else:
                now = time.time()
                self.prune_recent_submissions()
                self.update_queued_jobs()
                completions = self.check_for_completions()
                self.recent_completions.update(completions)
                self.prune_recent_completions()
                self._log(
                    "; ".join(
                        [
                            "{num_completions} jobs just completed",
                            "{num_pending} jobs in engine's queue",
                            "{num_recent_submissions} recent_submissions",
                        ]).format(
                            num_completions=len(completions),
                            num_pending=len(self.queued_jobs),
                            num_recent_submissions=len(self.recent_submissions),
                        )
                )
                if not self.cooling_down:
                    if self.tick_counter == 1:
                        self._log(
                            "First tick, clearing recent completions.")
                        self.recent_completions = {}
                    just_submitted_ids = self.submit_from_inbox()
                    self.recent_submissions.update(just_submitted_ids)
                    error_percent = self.get_error_percent()
                    lookback_start = datetime.datetime.now() - \
                            datetime.timedelta(seconds=self.error_lookback_time)
                    self._log(
                        "%s completions since %s, "
                        " error rate: %.1f" % (
                            len(self.recent_completions),
                            lookback_start,
                            error_percent
                        )
                    )
                    if error_percent > self.max_error_percent:
                        msg = "Exceeded max_error_percent of {}.".format(
                            self.max_error_percent)
                        self._log(msg)
                        if self.tick_counter == 1:
                            msg = ("First tick, ignoring high error rate.")
                        elif self.get_num_recent_errors() < \
                                self.error_lookback_thresh:
                            msg = (
                                "Below lookback error count threshold of"
                                " '{error_lookback_thresh}',"
                                " ignoring high error rate."
                            ).format(
                                error_lookback_thresh=self.error_lookback_thresh
                            )
                        else:
                            msg = "Pausing submissions for {} seconds.".format(
                                self.error_cooldown_time
                            )
                            self.cooling_down = True
                            self.error_cooldown_start = now
                        self._log(msg)
                else:
                    self._log("submission paused")
                    cooldown_finished = (now - self.error_cooldown_start) > \
                            self.error_cooldown_time
                    if cooldown_finished:
                        self.cooling_down = False
                        self.error_cooldown_start = None
                        self._log("resuming submissions")
                if self.do_cpu_gpu is True:
                    self.check_cpu_in_gpu(
                        'aspuru-samsung',
                        'aspuru-samsung-gpu',
                        'aspuru-samsung-gpu,aspuru-guzik,general,aspuru-samsung,serial_requeue',
                        'normal',
                        58,
                        3600)
                if self.max_error_count and \
                   (self.engine.get_error_count() > self.max_error_count):
                    print("STOPPING FROM TOO MANY ERRORS (details in log file)")
                    logging.fatal(self.engine.get_errors())
                    break

            self._log("{inbox_size} jobs waiting in inbox.".format(
                inbox_size=self.inbox_size))
            self._log(("Sleeping {period} seconds,"
                                 " next tick will happen at {next_tick_time}."
                                ).format(
                                    period=self.period,
                                    next_tick_time=(
                                        datetime.datetime.now() + \
                                        datetime.timedelta(seconds=self.period)
                                    )
                                ))
            time.sleep(self.period)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-r', '--root', dest='root', type=os.path.expanduser, default=".",
                        help="root path to make dirs")
    parser.add_argument('-e', '--engine', dest='engine', type=str, default="engines.SlurmEngine",
                        help="the submission engine module")
    parser.add_argument('-p', '--period', dest='period', type=int, default=DEFAULT_PERIOD_SECONDS,
                        help="number of seconds between updates (default: {0})".format(DEFAULT_PERIOD_SECONDS))
    parser.add_argument('--do_cpu_gpu', dest='do_cpu_gpu',  action='store_true', default=DEFAULT_DO_CPU_GPU,
                        help="whether to do cpu-gpu managing, default is {0}".format(DEFAULT_DO_CPU_GPU))
    parser.add_argument('-d', '--delay', dest='submit_delay', type=int, default=DEFAULT_SUBMIT_DELAY,
                        help="number of seconds between submitting jobs (default: {0})".format(DEFAULT_SUBMIT_DELAY))
    parser.add_argument('-M', '--maxjobs', dest='maxjobs', type=int, default=DEFAULT_MAX_QUEUE_LENGTH,
                        help="max total number of jobs (running or pending) (default={0})".format(DEFAULT_MAX_QUEUE_LENGTH))
    parser.add_argument('-m', '--maxpend', dest='maxpending', type=int, default=DEFAULT_MAX_PENDING_LENGTH,
                        help="max number of pending jobs (default={0})".format(DEFAULT_MAX_PENDING_LENGTH))
    help_str = "stop running are over N jobs with errors in queue (default N={0})".format(DEFAULT_MAX_ERROR_COUNT)
    parser.add_argument('-x', '--maxerrors', dest='maxerrors', type=int, default=DEFAULT_MAX_ERROR_COUNT,
                        help=help_str)
    parser.add_argument('-s', '--setup', action="store_true", default=False,
                        help="setup dirs in provided root dir (default=False)")
    parser.add_argument('-g', '--group_writable', action="store_true", default=False,
                        help="set all output files as group writable to match group of encompassing dir")
    parser.add_argument('-L', '--logfile', dest="logfilename", type=str,
                        default=DEFAULT_LOG_FILENAME,
                        help="log file name (default={0})".format(DEFAULT_LOG_FILENAME))
    parser.add_argument('-l', '--loglevel', dest='loglevel', type=int, default=DEFAULT_LOG_LEVEL,
                        help="log level (default: {0})".format(DEFAULT_LOG_LEVEL))
    parser.add_argument('--presubmit_checkpoint_filename',
                        help=("""The name of a 'presubmit' checkpoint file.
                              If this argument is specified, then a job dir in
                              inbox will not be submitted unless the file with
                              this name exists in the job dir. The intention of
                              this file is to ensure that job dirs are fully
                              copied before being submitted."""))
    parser.add_argument('--finished_checkpoint_filename',
                        help=("""The name of a 'finished' checkpoint file. If
                              specified, this file will be added to finished
                              job dirs to indicate that they are (1) finished,
                              and (2) copied to their final destination."""),
                        default=DEFAULT_FINISHED_CHECKPOINT_FILENAME)
    parser.add_argument('--global_lock_file',
                        type=os.path.expanduser,
                        help=("""Path to a 'global lock' file. If this file
                              exists then job_manager 'run' ticks will be
                              no-ops. The intended purpose of this lock
                              file is to allow workflow scripts to
                              temporarily lock the job_manager, to prevent
                              errors due to file copying."""))
    parser.add_argument('--engine_min_scan_interval', type=int,
                        help=("""Minimum time in seconds between scans of
                              an engine's queue. This is intended to prevent
                              overloading an engine's queue system with
                              queries."""),
                        default=DEFAULT_ENGINE_MIN_SCAN_INTERVAL)
    parser.add_argument('--submission_grace_period', type=int,
                        help=("""Time in seconds for which a just submitted
                              job will be considered as 'running'. The purpose
                              of this grace period is to allow for some slosh
                              between when a job is submitted, and when a 
                              job engine (e.g. slurm) reports a job as
                              'running'. If left empty, this will be set to
                              (2 * engine's minimum scan interval)""")),
    parser.add_argument('--max_error_percent', type=float,
                        help=("""Max percentage of failed jobs in the last
                              error_lookback seconds."""), default=50.0),
    parser.add_argument('--error_lookback_time', type=int,
                        help=("""Time in seconds to look back when calculating
                              error percent."""), default=(60 * 10)),
    parser.add_argument('--error_cooldown_time', type=int,
                        help=("""Time in seconds to delay before submitting more
                              jobs, when the error percentage has been
                              exceeded."""), default=1800),
    parser.add_argument('--error_lookback_thresh', type=int,
                        help=("""Minimum number of errors 
                              that can occur during error lookback window
                              before job_manager pauses for cooldown.
                              """), default=5),
    parser.add_argument('--scratch_dir',
                        type=os.path.expanduser,
                        help=("""Scratch dir for temp files.
                              If unset, will defautl to '<root>/.scratch'"""))
    clargs = parser.parse_args()

    logging.basicConfig(filename=os.path.join(clargs.root, clargs.logfilename),
                        level=clargs.loglevel,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%m-%d %H:%M')

    engine_module_name, engine_class_name = clargs.engine.split(".")
    engine_module = importlib.import_module(engine_module_name)
    engine = getattr(engine_module, engine_class_name)(
        min_scan_interval=clargs.engine_min_scan_interval)

    if engine_class_name == 'LocalEngine':
        submit_delay = DEFAULT_SUBMIT_DELAY_LOCAL
    else:
        submit_delay = clargs.submit_delay

    submission_grace_period = clargs.submission_grace_period
    if not submission_grace_period:
        submission_grace_period = 2 * clargs.engine_min_scan_interval

    scratch_dir = clargs.scratch_dir or os.path.join(clargs.root, '.scratch')

    jm = JobManager(
        clargs.root,
        engine,
        period=clargs.period,
        submit_delay=submit_delay,
        max_queue_len=clargs.maxjobs,
        max_pending_len=clargs.maxpending,
        max_error_count=clargs.maxerrors,
        setup=clargs.setup,
        do_cpu_gpu=clargs.do_cpu_gpu,
        output_group_writable=clargs.group_writable,
        presubmit_checkpoint_filename=clargs.presubmit_checkpoint_filename,
        finished_checkpoint_filename=clargs.finished_checkpoint_filename,
        submission_grace_period=submission_grace_period,
        global_lock_file=clargs.global_lock_file,
        max_error_percent=clargs.max_error_percent,
        error_lookback_time=clargs.error_lookback_time,
        error_cooldown_time=clargs.error_cooldown_time,
        error_lookback_thresh=clargs.error_lookback_thresh,
        scratch_dir=scratch_dir,
    )
    if not clargs.setup:
        jm._log("Starting Job Manager: use ctrl-C to quit")
        jm._log("Watching inbox: " + jm.inbox_dir)
        jm.run()
