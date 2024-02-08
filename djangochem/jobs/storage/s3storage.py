import datetime
import functools
import json
import os
import re
import shutil
import signal
import time
from distutils.dir_util import copy_tree

import boto3
import botocore


from .signalhandler import SignalHandler

TIMESTAMP_SCHEME = "%y-%m-%d-%H:%M"
JOB_NAME_SCHEME = "{key}_{timestamp}_{name}"
DEFAULT_BUCKET_NAME = 'calculario.default'


class S3WriteDir(object):

    def __init__(self, key, name, priority=5, inbox_job_dir='inbox', bucket_name=DEFAULT_BUCKET_NAME):
        s3 = boto3.resource('s3')
        self.bucket = s3.Bucket(bucket_name)
        # throws error if bucket doesn't exist
        try:
            s3.meta.client.head_bucket(Bucket=bucket_name)
        except botocore.exceptions.ClientError:
            s3.create_bucket(Bucket=bucket_name)
            self.bucket = s3.Bucket(bucket_name)
        self.name = name
        now_str = datetime.datetime.strftime(datetime.datetime.utcnow(), TIMESTAMP_SCHEME)
        job_name = JOB_NAME_SCHEME.format(timestamp=now_str,
                                          key=key,
                                          name=name)
        self.job_path = os.path.join(inbox_job_dir, job_name)
        if list(self.bucket.objects.filter(Prefix=self.job_path)):
            raise Exception("S3 Key exists at {}".format(self.job_path))

    def dump_json(self, filename, info):
        json_str = json.dumps(info)
        path = os.path.join(self.job_path, filename)
        self.bucket.Object(path).put(Body=json_str)
        return path

    def dump_file(self, filename, contents):
        path = os.path.join(self.job_path, filename)
        self.bucket.Object(path).put(Body=contents)
        return path

    def make_subdir(self, dirname):
        """ added for consistency with file storage, but dirs don't exist in s3"""
        pass


    def __repr__(self):
        return "<S3WriteDir " + self.job_path + ">"


class S3ReadDir(S3WriteDir):

    def __init__(self, job_dir_path, name='', error_path="error", archive_path="archive", bucket_name=DEFAULT_BUCKET_NAME):
        s3 = boto3.resource('s3')
        self.bucket = s3.Bucket(bucket_name)
        # throws error if bucket doesn't exist
        s3.meta.client.head_bucket(Bucket=bucket_name)
        self.name = name
        self.job_path = job_dir_path
        if not list(self.bucket.objects.filter(Prefix=self.job_path)):
            raise Exception("No S3 Keys exists with prefix {}".format(self.job_path))

        self.error_path = error_path
        self.archive_path = archive_path

    def load_json(self, filename):
        json_filename = os.path.join(self.job_path, filename)
        json_str = self.bucket.Object(json_filename).get()['Body'].read()
        info = json.loads(json_str.decode('utf8'))
        return info

    def copy_to_path(self, dest_dir):
        for objectsummary in self.bucket.objects.filter(Prefix=self.job_path):
            slash_and_filename = objectsummary.key.split(self.job_path)[1]
            dest_path = dest_dir + slash_and_filename
            # ensure any intermediate sub-directories get built
            sub_dir = os.path.dirname(dest_path)
            if sub_dir != dest_dir:
                os.makedirs(sub_dir, exist_ok=True)
            objectsummary.Object().download_file(dest_path)

    def file_exists(self, filename):
        objs = list(self.bucket.objects.filter(Prefix=filename))
        return len(objs) > 0 and objs[0].key == filename

    def dir_exists(self, dirname):
        """ no such thing as a dir, so this means, if there are keys with the prefix"""
        if not dirname.endswith('/'):
            dirname += '/'
        objs = list(self.bucket.objects.filter(Prefix=dirname))
        return len(objs) > 0 and objs[0].key.startswith(dirname) and objs[0].key != dirname

    def list_dir(self, dir_only=False):
        out = []
        for o in self.bucket.objects.filter(Prefix=self.job_path):
            if not dir_only or self.dir_exists(o.key):
                out.append(o.key)
        return out

    def mark_error(self):
        return self._move_to(self.error_path)

    def archive(self):
        return self._move_to(self.archive_path)

    def _move_to(self, dest_dir):
        job_dir_name = os.path.basename(self.job_path)
        clean_dest_dir = os.path.join(dest_dir, job_dir_name)
        dest_dir = clean_dest_dir
        if self.dir_exists(dest_dir):
            i = 2
            while self.dir_exists(dest_dir):
                dest_dir = clean_dest_dir + "_{0}".format(i)
                i += 1
            self.logger.warning("Job dir already exists in archive dir, renaming to {0}".format(dest_dir))
        with SignalHandler([signal.SIGINT]) as handler:
            for objectsummary in self.bucket.objects.filter(Prefix=self.job_path):
                slash_and_filename = objectsummary.key.split(self.job_path)[1]
                dest_path = dest_dir + slash_and_filename
                # ensure any intermediate sub-directories get built
                bucket_and_src_path = os.path.join(self.bucket.name, objectsummary.key)
                self.bucket.Object(dest_path).copy_from(CopySource=bucket_and_src_path)
                objectsummary.Object().delete()
            if handler.interrupted:
                msg = "Terminated with signal {0}".format(handler.sig_used)
                raise(Exception(msg))

    def __repr__(self):
        return "<S3ReadDir " + self.job_path + ">"
