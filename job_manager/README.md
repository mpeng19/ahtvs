/shareJob Manager

Author: Tim Hirzel <hirzel@chemistry.harvard.edu>
Date: April 2014

Intro
-----

job_manager.py is an app that runs on odyssey and allows you to send jobs and get answers via files.  It is a process that watches an inbox directory and calls sbatch on job files that it finds.  It maintains three directories: inbox, pending, and complete.  It will look for folders inside inbox that contain a file job.sh inside of them.  For example `inbox/job_10/job.sh`.

Jobs will be pulled from the inbox in alphabetical order.  Use this feature to prioritize your jobs.  An example system might be to simple number your folders where lower numbers are higher priority, such as `inbox/100/job.sh`, `inbox/110/job.sh`, and `inbox/140/job.sh`.

Any files that you put in your job folders will travel with that folder from the inbox, into the pending and finally into the completed dir.  The job manager will write some helper files into this directory as well.

In the case that your job reaches an error state. It will end up in the completed folder with a special job_error.txt file.  This file will have the error line from the queue written into it.  You will probably have slurm output with more error info as well, but there are cases where errors occur in the queue before slurm can even get to it.

The job_manager has no state.  If your job_manager is killed, don't worry.  Just restart it and your directories should continue to progress through the directories.  The only thing you might lose if it was down for more than a day and you had a long running jobs that ended in an error.  The job manager may miss detecting an error from the queue if it was purged before reading the queue.

To setup the job_manager
-------------------------

0. get this job_manager folder onto your Odyssey drive.  either with git clone or just scp the dir.

10. create a new directory somewhere in your odyssey share dir that is accessible on both machines
 This will become the "root" of the all activity with the inbox, pending, and complete subfolers in it
 for this tutorial, let's assume you do:

        odyssey$ mkdir ~/jobs

    This directory should be one that you can see from your local machine.  You will be putting job directories into the inbox directory and getting completed jobs from the completed directory.  We recommend that you use ssfs on your local machine to mount this folder locally to get the most enjoyment of the job manager.  You can use sftp or scp to move job folders in and out as well if you prefer.

20. Load the python2.7.3 module

        odyssey$ module load hpc/python-2.7.3

30. cd into the job_manager dir.  (this is optional, but if you don't, you must compensate by providing the full path to the job_manager.py file in the future steps)

        odyssey$ cd job_manager

35. familiarize yourself with the job_manager options and defalts by looking at the help

        odyssey$ python job_manager.py --help

4. setup the share directory (create the inbox, pending, and complete directories)

        odyssey$ python job_manager.py --root=~/jobs --setup

6. start it up!

        odyssey$ python job_manager.py --root=~/jobs

    you will see a startup message like this:

        [hirzel@rclogin12 ~]$ python job_manager/job_manager.py --root=share
        Starting Job Manager: use ctrl-C to quit
        Watching inbox: /n/home13/hirzel/share/inbox
        ....

    The dots will print out in rows of 20 as a heartbeat for the app.  One dot prints every time
    the app checks up on the jobs.  This occurs every number of seconds specified to the job_manager with the --period option (default 60).  It is safe to use Ctrl-C to stop the app.  It has signal handler to prevent it stopping in a bad spot.

7. put some directories containing job.sh files into the inbox.  enjoy!

Running the Job Manager locally
--------------------------------

In this mode, your job.sh script will be executed simply as a bash script with "bash job.sh"
any stderr or stdout will be directed to `job_output.txt`

To do this set command line option  `-e` to `engines.LocalEngine`
for example:
        python job_manager.py -r ~/local_jobs -s # setup local job dirs
        python job_manager.py -r ~/local_jobs -e engines.LocalEngine -p2
