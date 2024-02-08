import sys
import os
import shutil
import tarfile
import random
import django
import time
from django.core.management import call_command


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

jobbuild_dir = os.getenv("PIPE_BUILD_DIR", "build")
jobtemp_dir = os.getenv("PIPE_TMP_DIR", "tmp")


# slurm_job_id = os.getenv("SLURM_JOB_ID", None)
# if (slurm_job_id is None) or (slurm_job_id == ''):
#     slurm_job_id = str(random.randint(1000000, 9999999)) # '360475'#

# addsmiles sf_lib ~/data/reaxys.smi -t reaxys4 --settings=djangochem.settings.denn_reaxys4_rc

def _move_dir_to_tar_local(source_dir=None, dest_tar=None, compression='gz'):
    with tarfile.open(
            name=dest_tar,
            mode='w:{}'.format(compression)
    ) as tar_handle:
        tar_handle.add(source_dir, arcname='')
    shutil.rmtree(source_dir)


def buildjob(chemconf, batch):
    # os.mknod(os.path.join(jobbuild_dir, LOCK_FILE_NAME))
    # jobtemp_temp_dir = os.path.join(jobtemp_dir, str(time.time()))
    jobtemp_temp_dir = jobtemp_dir
    if not os.path.exists(jobtemp_temp_dir):
        os.makedirs(jobtemp_temp_dir)
    if batch>1:
        call_command('buildjobs', project_name, jobtemp_temp_dir, config=chemconf, batchsize=batch)
    else:
        call_command('buildjobs', project_name, jobtemp_temp_dir, config=chemconf)

    for f in os.scandir(jobtemp_temp_dir):
        if f.is_dir():
            tar_name = "{}.tar.gz".format(str(f.name))
            tar_name_temp = tar_name + '.temp'
            tmp_tar = os.path.join(jobtemp_temp_dir, tar_name_temp)
            _move_dir_to_tar_local(f.path, tmp_tar, 'gz')
            tar_inbox_temp_path = os.path.join(jobbuild_dir, tar_name_temp)
            shutil.move(tmp_tar, tar_inbox_temp_path)
            os.rename(tar_inbox_temp_path, os.path.join(jobbuild_dir, tar_name))

    shutil.rmtree(jobtemp_temp_dir)
    # os.remove(os.path.join(jobbuild_dir, LOCK_FILE_NAME))

for i in range(ITERATIONS):
    print("Iteration {} of {} started...".format(i,ITERATIONS))

    # Parse jobs moved to separate script
    #extract_parse_completed()

    call_command('requestjobs', project_name, 'conformer_clean_dftb', tag=tag, limit=100)
    buildjob('conformer_clean_dftb',20)

    call_command('requestjobs', project_name, 'dftb_d3_opt', parent_config='conformer_clean_dftb', tag=tag, limit=100)
    buildjob('dftb_d3_opt',100)

    call_command('requestjobs', project_name, 'dftb_tddft', parent_config='dftb_d3_opt', tag=tag, limit=100)
    buildjob('dftb_tddft',100)

    print("Sleeping...")
    time.sleep(PERIOD)
