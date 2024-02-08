from distutils.dir_util import copy_tree
import datetime
import functools
import json
import logging
import os
import re
import shutil
import signal
import time

from .signalhandler import SignalHandler


MAX_DELETE_TRIES = 10
JOB_NAME_SEPARATOR = "_"  # must match below
ARCHIVE_DIR_NAME_SCHEME = "%y-%m-%d"
JOB_NAME_SCHEME = "{num:04d}_{name}_{key}_{timestamp}"


class WriteDir(object):

    def __init__(self, key, name, priority=5, inbox_job_dir='', number_hint=None):
        self.job_path = None
        self.mynumber = None
        while self.mynumber is None or os.path.exists(self.job_path):
            self.job_path = self._make_dirpath(key, name, priority, inbox_job_dir)
        os.mkdir(self.job_path)

    def _make_dirpath(self, key, name, priority, inbox_job_dir):
        if self.mynumber is None:
            self.mynumber = self._get_next_job_number(inbox_job_dir, priority)
        else:
            self.mynumber += 1
        self.name = name
        job_name = JOB_NAME_SCHEME.format(num=self.mynumber,
                                          key=key,
                                          name=name,
                                          timestamp=int(time.time()))
        return os.path.join(inbox_job_dir, job_name)


    def _get_next_job_number(self, inbox_job_dir, priority):
        """
        scan inbox and create a new job number that is one greater than the max present
        empty dir or non-numbered dirs will lead to return val of 0
        """
        parsenum = re.compile("\d+")
        def num_cmp(a, b):
            try:
                anum = int(parsenum.findall(a)[0])
            except:
                return -1
            try:
                bnum = int(parsenum.findall(b)[0])
            except:
                return 1
            return anum > bnum

        dirs = [d for d in os.listdir(inbox_job_dir) if d.startswith(str(priority))]
        dirs.sort(key=functools.cmp_to_key(num_cmp))
        try:
            last = dirs[-1]
            curmax = int(parsenum.findall(last)[0])
            return curmax + 1
        except IndexError:
            return priority * 1000

    def dump_json(self, filename, info):
        with open(os.path.join(self.job_path, filename), 'w') as fd:
            json.dump(info, fd)

    def dump_file(self, filename, contents):
        full_file_path = os.path.join(self.job_path, filename)
        file_dir = os.path.dirname(full_file_path)
        if not os.path.exists(file_dir):
            os.makedirs(file_dir)
        with open(full_file_path, 'w') as fd:
            fd.write(contents)

    def make_subdir(self, dirname):
        os.mkdir(os.path.join(self.job_path, dirname))

    def __repr__(self):
        return "<WriteDir " + self.job_path + ">"


class ReadDir(WriteDir):

    def __init__(self, job_dir_path, name='', error_path=None, archive_path=None):
        if not os.path.isdir(job_dir_path):
            raise Exception("{} is not a directory".format(job_dir_path))
        self.job_path = job_dir_path
        self.name = name
        self.error_path = error_path
        self.archive_path = archive_path

    def list_dir(self, dir_only=False):
        """return contents of job dir as full paths"""
        out = []
        for p in os.listdir(self.job_path):
            full_path = os.path.join(self.job_path, p)
            if not dir_only or os.path.isdir(full_path):
                out.append(full_path)
        return out

    def load_json(self, filename):
        info_filename = os.path.join(self.job_path, filename)
        with open(info_filename, 'r') as f:
            info = json.load(f)
        return info

    def copy_to_path(self, dest_path):
        copy_tree(self.job_path, dest_path)

    def file_exists(self, filename):
        return os.path.isfile(filename)

    def dir_exists(self, dirname):
        return os.path.isdir(dirname)

    def mark_error(self):
        if self.error_path is not None:
            return self._move_to(self.job_path, self.error_path)

    def archive(self):
        if self.archive_path is not None:
            return self._move_to(self.job_path, self.archive_path, use_clock_subdir=True)

    def _move_to(self, src_path, dest_dir, use_clock_subdir=False):
        cur_path = src_path
        base_src_path = os.path.basename(src_path)
        if use_clock_subdir:
            now_str = datetime.datetime.strftime(datetime.datetime.utcnow(), ARCHIVE_DIR_NAME_SCHEME)
            current_day_worker_dir = self.name + "_" + now_str
            base_dest_path = os.path.join(dest_dir, current_day_worker_dir, base_src_path)
        else:
            base_dest_path = os.path.join(dest_dir, base_src_path)
        dest_path = base_dest_path
        i = 1
        if os.path.exists(dest_path):
            while os.path.exists(dest_path):
                dest_path = base_dest_path + "_DUP{0}".format(i)
                i += 1
            logging.warning("Job dir already exists in archive dir, renaming to {0}".format(dest_path))
        with SignalHandler([signal.SIGINT]) as handler:
            try:
                if os.path.isdir(cur_path):
                    shutil.copytree(cur_path, dest_path)
                else:
                    # check that the dest dir is actually there in the case of moving file
                    dest_dir = os.path.dirname(dest_path)
                    if not os.path.exists(dest_dir):
                        os.makedirs(dest_dir)
                    shutil.copyfile(cur_path, dest_path)
                deleted = False
                attempts = 0
                while not deleted:
                    try:
                        if os.path.isdir(cur_path):
                            shutil.rmtree(cur_path)
                        else:
                            os.remove(cur_path)
                        deleted = True
                    except (IOError, OSError):
                        attempts += 1
                        time.sleep(0.01)
                        if attempts >= MAX_DELETE_TRIES:
                            raise
            except (IOError, OSError):
                msg = "Shutil error when moving the file"
                logging.exception(msg)
            if handler.interrupted:
                msg = "Terminated with signal {0}".format(handler.sig_used)
                logging.warning(msg)
                raise Exception(msg)

    def __repr__(self):
        return "<ReadDir " + self.job_path + ">"
