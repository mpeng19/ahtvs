import unittest
import os
import shutil

from job_manager import job_manager
from job_manager.engines import Engine

MYDIR = os.path.dirname(__file__)
TEST_SHARE_DIR = "test_share"


class TestEngine(Engine):
    COMPLETED_MARKER = "completed"

    def __init__(self):
        self.jobs = []

    def submit_job(self, path, name="noname"):
        self.jobs.append(["user", "1", name, "completed"])
        return "1"

    def get_pending(self):
        return self.jobs

    def get_errors(self):
        return [("user", "1", "blah", "test error")]

    def get_all(self):
        return self.jobs


class TestJobManager(unittest.TestCase):

    def setUp(self):
        self.test_engine = TestEngine()
        self.testdir = os.path.join(MYDIR, TEST_SHARE_DIR)
        os.makedirs(self.testdir)
        self.jm = job_manager.JobManager(self.testdir,
                                         self.test_engine,
                                         submit_delay=0,
                                         max_queue_len=1,
                                         setup=True)

    def make_job(self, job_id):
        test_job_dirpath = os.path.join(self.testdir,
                                        job_manager.INBOX_NAME,
                                        job_id
                                        )
        os.makedirs(test_job_dirpath)
        test_job_submit_path = os.path.join(test_job_dirpath,
                                            job_manager.DEFAULT_JOB_FILENAME)
        f = open(test_job_submit_path, 'w')
        f.write("test job")

    def test_submit(self):
        job_id = "0-job_id"
        self.make_job(job_id)
        self.jm.check_inbox()
        self.assertEquals(len(self.test_engine.jobs), 1)
        test_pending_dirpath = os.path.join(self.testdir,
                                            job_manager.PENDING_NAME)
        self.assertIn(job_id, os.listdir(test_pending_dirpath))
        return job_id

    def test_complete(self):
        job_id = "0-job_id"
        self.make_job(job_id)
        self.jm.check_inbox()
        self.jm.check_for_completed()
        self.jm.check_for_completed()
        test_complete_dirpath = os.path.join(self.testdir,
                                             job_manager.COMPLETED_NAME)
        self.assertIn(job_id, os.listdir(test_complete_dirpath))

    def test_submit_twice_hit_limit(self):
        self.make_job("2")
        self.make_job("1")
        self.jm.check_inbox()
        self.jm.check_inbox()
        test_pending_dirpath = os.path.join(self.testdir,
                                            job_manager.PENDING_NAME)
        self.assertIn("1", os.listdir(test_pending_dirpath))
        test_inbox_dirpath = os.path.join(self.testdir,
                                          job_manager.INBOX_NAME)
        self.assertIn("2", os.listdir(test_inbox_dirpath))

    def test_submit_twice_then_release(self):
        self.test_submit_twice_hit_limit()
        self.jm.check_for_completed()
        self.jm.check_for_completed()
        self.test_engine.jobs = []  # clear jobs as if done
        self.jm.check_inbox()
        test_pending_dirpath = os.path.join(self.testdir,
                                            job_manager.PENDING_NAME)
        self.assertIn("2", os.listdir(test_pending_dirpath))
        test_complete_dirpath = os.path.join(self.testdir,
                                             job_manager.COMPLETED_NAME)
        self.assertIn("1", os.listdir(test_complete_dirpath))

    def test_error(self):
        self.make_job("1-err")
        self.jm.check_inbox()

        err_count = self.jm.check_for_errors()
        self.assertEquals(err_count, 1)
        test_complete_dirpath = os.path.join(self.testdir,
                                             job_manager.COMPLETED_NAME)
        self.assertIn("1-err", os.listdir(test_complete_dirpath))
        self.assertIn(job_manager.DEFAULT_ERROR_FILENAME,
                      os.listdir(test_complete_dirpath + "/1-err"))
        with open(test_complete_dirpath + "/1-err/" + job_manager.DEFAULT_ERROR_FILENAME, 'r') as f:
            print(f.read())

    def test_pid_file(self):
        self.make_job("priority1")
        self.jm.check_inbox()
        self.assertEquals(len(self.test_engine.jobs), 1)
        test_pending_dirpath = os.path.join(self.testdir,
                                            job_manager.PENDING_NAME)
        self.assertIn(job_manager.DEFAULT_JOB_ID_FILENAME,
                      os.listdir(test_pending_dirpath+"/priority1"))

    def tearDown(self):
        shutil.rmtree(self.testdir)

if __name__ == "__main__":
    unittest.main()
