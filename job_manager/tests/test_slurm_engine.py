import unittest
import os
import shutil
import time
from job_manager.engines import SlurmEngine
import random

test_script = """#!/bin/bash
#SBATCH -p aspuru-guzik # partition
sleep 2
touch complete.txt
"""

JOB_DIR = "TEST_JOB_0"
JOB_NAME = "job.sh"

MAX_WAIT = 120


class TestSlurmEngine(unittest.TestCase):

    def setUp(self):
        self.job_dir = JOB_DIR+str(random.randint(0, 1000))
        os.mkdir(self.job_dir)
        self.job_path = os.path.join(self.job_dir, JOB_NAME)
        self.engine = SlurmEngine()
        with open(self.job_path, 'w') as script:
            script.write(test_script)

    def test_job(self):
        print(self.engine.submit_job(self.job_path))
        d = 0
        while "complete.txt" not in os.listdir(self.job_dir):
            time.sleep(5)
            d += 5
            if d > MAX_WAIT:
                raise Exception("time out waiting")

    def test_queue_len(self):
        print(self.engine.submit_job(self.job_path))
        d = 0
        while self.engine.get_pending_count() != 1:
            time.sleep(5)
            d += 5
            if d > MAX_WAIT:
                raise Exception("time out waiting")

    def test_error_count(self):
        self.assertEquals(self.engine.get_error_count(), 0)

    def tearDown(self):
        shutil.rmtree(self.job_dir)

if __name__ == "__main__":
    unittest.main()
