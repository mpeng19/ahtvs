import glob
import os
import shutil
import tarfile
import django
import time
from django.core.management import call_command
import random


# setup the django settings file.  Change this to use the settings file that connects you to your desired database
# os.environ["DJANGO_SETTINGS_MODULE"] = sys.argv[1]
# this must be run to setup access to the django settings and make database access work etc.
django.setup()

ITERATIONS = int(os.getenv("PIPE_BUILD_ITER", 10000))
PERIOD = 1
#LIMIT = 200
#BATCH = 200
tag = [os.getenv("PIPE_TAG", "d1")]
project_name = os.getenv("PIPE_PROJECT", "demo")

jobcomplete_dir = os.getenv("PIPE_COMPL_DIR", "build")
jobtemp_dir = os.getenv("PIPE_TMP_DIR", "tmp")
jobarchive_dir = os.getenv("PIPE_ARCH_DIR", "tmp")


def extract_parse_completed():
    jobs = glob.glob(os.path.join(jobcomplete_dir, '*.tar.gz'), recursive=False)

    if len(jobs)>0:
        job=random.choice(jobs)
        #for job in jobs[:1]:
        try:
            shutil.move(job, jobtemp_dir)
        except:
            print("Conflict!")
            return

        jobs_tars_local = glob.glob(os.path.join(jobtemp_dir, '*.tar.gz'), recursive=False)
        for job_tar in jobs_tars_local:
            job_dir = os.path.join(jobtemp_dir, job_tar.rstrip('.tar.gz'))
            with tarfile.open(name=job_tar,
                              mode='r:gz') as handle:
                handle.extractall(path=job_dir)

            if os.path.exists(job_tar):
                os.remove(job_tar)

            call_command('parsejobs', project_name, jobtemp_dir, root_path=jobarchive_dir)

            if os.path.exists(job_dir):
                shutil.rmtree(job_dir)
    else:
        print("WARNING: no jobs to parse")

if not os.path.exists(jobtemp_dir):
    os.makedirs(jobtemp_dir)

    extract_parse_completed()

#cleanup
if os.path.exists(jobtemp_dir):
    shutil.rmtree(jobtemp_dir)