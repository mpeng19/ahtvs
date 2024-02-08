import os
import unittest
import tempfile

try:
    import boto3
    BOTO_AVAIL = True
    from ..storage.s3storage import S3WriteDir, S3ReadDir
except ImportError:
    BOTO_AVAIL = False

TEST_BUCKET = "calculario.unittest"
JSON_NAME = "testfile.json"
FILE_NAME = "test.txt"

TEST_DICT = {"hello": "world"}
TEST_STRING = "hello there world!"

@unittest.skipIf(not BOTO_AVAIL, reason="requires boto3")
class TestS3Storage(unittest.TestCase):
    def setUp(self):
        s3 = boto3.resource('s3')
        self.bucket = s3.Bucket(TEST_BUCKET)
        self.writer = S3WriteDir("KEY", "test", "inbox", bucket_name=TEST_BUCKET)

    def _add_file(self, filename=FILE_NAME, content=TEST_STRING):
        self.writer.dump_file(filename, content)

    def _add_json_file(self):
        self.writer.dump_json(JSON_NAME, TEST_DICT)

    def test_write(self):
        self._add_json_file()
        self.assertEquals(len(list(self.bucket.objects.all())), 1)

    def test_read(self):
        self._add_json_file()
        reader = S3ReadDir(self.writer.job_path, bucket_name=TEST_BUCKET)
        self.assertEquals(TEST_DICT, reader.load_json(JSON_NAME))

    def test_dir_dump(self):
        files = [('a.txt', 'aaa'),
                 ('b.txt', 'bbb'),
                 ('c.txt', 'ccc'),
                 ('subdir/d.txt', 'ddd')]
        for name, content in files:
            self._add_file(name, content)
        reader = S3ReadDir(self.writer.job_path, bucket_name=TEST_BUCKET)

        subdir = os.path.join(self.writer.job_path, "subdir")
        self.assertTrue(reader.dir_exists(subdir))

        with tempfile.TemporaryDirectory() as tempdir:
            reader.copy_to_path(tempdir)
            self.assertTrue(os.path.isdir(os.path.join(tempdir, 'subdir')))
            for name, content in files:
                temp_path = os.path.join(tempdir, name)
                self.assertEquals(content, open(temp_path).read())

    def test_error_move(self):
        files = [('a.txt', 'aaa'),
                 ('b.txt', 'bbb'),
                 ('c.txt', 'ccc'),
                 ('subdir/d.txt', 'ddd')]
        for name, content in files:
            self._add_file(name, content)
        reader = S3ReadDir(self.writer.job_path, bucket_name=TEST_BUCKET)
        reader.mark_error()
        for name, content in files:
            old_path = os.path.join(self.writer.job_path, name)
            new_path = os.path.join(*["error"] + old_path.split("/")[1:])
            self.assertFalse(reader.file_exists(old_path))
            self.assertTrue(reader.file_exists(new_path))

    def test_archive_move(self):
        files = [('a.txt', 'aaa'),
                 ('b.txt', 'bbb'),
                 ('c.txt', 'ccc'),
                 ('subdir/d.txt', 'ddd')]
        for name, content in files:
            self._add_file(name, content)
        reader = S3ReadDir(self.writer.job_path, bucket_name=TEST_BUCKET)
        reader.archive()
        for name, content in files:
            old_path = os.path.join(self.writer.job_path, name)
            new_path = os.path.join(*["archive"] + old_path.split("/")[1:])
            self.assertFalse(reader.file_exists(old_path))
            self.assertTrue(reader.file_exists(new_path))

    def tearDown(self):
        for key in self.bucket.objects.all():
            key.delete()
