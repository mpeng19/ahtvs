import os

import django
import sys
# setup the djhango settings file.  Change this to use the settings file that connects you to your desired database
import shutil

import subprocess

import time

# <<<<<<< HEAD:pipeline_template/pipeline_template.py
os.environ["DJANGO_SETTINGS_MODULE"] = "djangochem.settings.postgres_docker"
# =======
# CHANGE ME
os.environ["DJANGO_SETTINGS_MODULE"] = "djangochem.settings.docker"
# >>>>>>> b282143d39b9cac5ebcb30743332cd16b2962ffd:examples/pipeline_template/pipeline_template.py
# this must be run to setup access to the django settings and make database access work etc.
django.setup()

from django.core.management import call_command

LOCK_FILE_NAME = "job_lock_file"
PERIOD = 500
LIMIT = 20
PM7_LIMIT = 1000
jobbuild_dir = "/Users/danieltabor/job_storage/gate_molecules/gate_molecules_jobdirs_build"
jobcomplete_dir = "/Users/danieltabor/job_storage/gate_molecules/gate_molecules_jobdirs_completed"
jobarchive_dir = "/Users/danieltabor/job_storage/gate_molecules/gate_molecules_jobdirs_parsed"
project_name = 'postgres_docker'

#specified_tag = ["symmetric_azulene_gen1"]

iterations = 0

def buildjob(chemconf):
    call_command('buildjobs', project_name, jobbuild_dir, config=chemconf)
    subprocess.call("mv /Users/danieltabor/job_storage/gate_molecules/gate_molecules_jobdirs_build/job_lock_file /Users/danieltabor/ody_gate_mols_jobs/inbox/", shell=True)
    subprocess.call("mv /Users/danieltabor/job_storage/gate_molecules/gate_molecules_jobdirs_build/* /Users/danieltabor/ody_gate_mols_jobs/inbox/", shell=True)
    subprocess.call("mv /Users/danieltabor/ody_gate_mols_jobs/inbox/job_lock_file /Users/danieltabor/job_storage/gate_molecules/gate_molecules_jobdirs_build/job_lock_file", shell=True)


def buildjob_batch(chemconf,batch_size):
    call_command('buildjobs', project_name, jobbuild_dir, config=chemconf, batchsize=batch_size)
    subprocess.call("mv /Users/danieltabor/job_storage/gate_molecules/gate_molecules_jobdirs_build/job_lock_file /Users/danieltabor/ody_gate_mols_jobs/inbox/", shell=True)
    subprocess.call("mv /Users/danieltabor/job_storage/gate_molecules/gate_molecules_jobdirs_build/* /Users/danieltabor/ody_gate_mols_jobs/inbox/", shell=True)
    subprocess.call("mv /Users/danieltabor/ody_gate_mols_jobs/inbox/job_lock_file /Users/danieltabor/job_storage/gate_molecules/gate_molecules_jobdirs_build/job_lock_file", shell=True)

while True:
    print("Working...")
    iterations = iterations + 1
    print("Iteration number:",iterations)

    subprocess.call("mv /Users/danieltabor/ody_gate_mols_jobs/completed/* /Users/danieltabor/job_storage/gate_molecules/gate_molecules_jobdirs_completed/", shell=True)

    #Parse jobs
    call_command('parsejobs', project_name, jobcomplete_dir, root_path=jobarchive_dir)

    call_command('requestjobs', project_name, 'conformer', limit=5)
    buildjob('conformer')

    call_command('requestjobs', project_name, 'b3lyp_6-31gs_opt_qchem', parent_config='conformer', limit=LIMIT,tag=specified_tag)
    buildjob('b3lyp_6-31gs_opt_qchem')

    call_command('requestjobs', project_name, 'b3lyp_6-31gs_sp_singlet_vertical_excited_qchem', parent_config='b3lyp_6-31gs_opt_qchem', limit=LIMIT)
    buildjob('b3lyp_6-31gs_sp_singlet_vertical_excited_qchem')

    call_command('requestjobs', project_name, 'b3lyp_6-31gs_sp_triplet_vertical_excited_qchem', parent_config='b3lyp_6-31gs_opt_qchem', limit=LIMIT)
    buildjob('b3lyp_6-31gs_sp_triplet_vertical_excited_qchem')

    call_command('requestjobs', project_name, 'b3lyp_6-31gs_opt_s1', parent_config='b3lyp_6-31gs_opt_qchem', limit=LIMIT)
    buildjob('b3lyp_6-31gs_opt_s1')

    call_command('requestjobs', project_name, 'b3lyp_6-31gs_opt_s2', parent_config='b3lyp_6-31gs_opt_qchem', limit=LIMIT)
    buildjob('b3lyp_6-31gs_opt_s2')

    call_command('requestjobs', project_name, 'b3lyp_6-31gs_opt_t1_qchem', parent_config='b3lyp_6-31gs_opt_qchem', limit=LIMIT)
    buildjob('b3lyp_6-31gs_opt_t1_qchem')


    print("Sleeping...")
    time.sleep(PERIOD)
