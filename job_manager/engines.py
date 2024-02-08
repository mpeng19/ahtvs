import subprocess
import os
import getpass
import logging
import time

SLURM_RUN_STATES = "CONFIGURING,COMPLETING,PENDING,RUNNING,RESIZING,SUSPENDED"
SLURM_ERROR_STATES = "CANCELLED,FAILED,TIMEOUT,PREEMPTED,NODE_FAIL"
SLURM_COMPLETED_STATES = "COMPLETED"

TORQUE_RUN_STATES = "E,H,Q,R,T,W,S"
TORQUE_ERROR_STATES =  "Error"
TORQUE_COMPLETED_STATES = "C"

SQUEUE_CMD = 'squeue  --noheader --format="%u %i %j %S %P %q %T"'
SACCT_CMD  = 'sacct --noheader  --format=user%30,jobid%12,jobname%100,start,partition%100,qos%20,state'

QSTAT_CMD = 'qstat -u bombarel'


MAX_ERRORS = 100

# from http://stackoverflow.com/questions/4814970/subprocess-check-output-doesnt-seem-to-exist-python-2-6-5
if "check_output" not in dir(subprocess):  # duck punch it in!
    def f(*popenargs, **kwargs):
        if 'stdout' in kwargs:
            raise ValueError('stdout argument not allowed, it will be overridden.')
        process = subprocess.Popen(stdout=subprocess.PIPE, *popenargs, **kwargs)
        output, unused_err = process.communicate()
        retcode = process.poll()
        if retcode:
            cmd = kwargs.get("args")
            if cmd is None:
                cmd = popenargs[0]
            raise subprocess.CalledProcessError(retcode, cmd)
        return output
    subprocess.check_output = f


class Engine:
    # if your engine supports pending only, set this to the string
    # that specifies a job is pending in the return of get_pending
    # this state would be the last item in any row of the list returned
    # from get_pending
    PENDING_MARKER = None

    def __init__(self, *args, min_scan_interval=None, **kwargs):
        self.min_scan_interval = min_scan_interval

    def submit_job(self, path, name=None):
        raise NotImplementedError

    def get_pending(self):
        """
        jobs pending completion as list of lists, where sublists have ordering such as
        the first two items are [user,jobid,...] and last column is state
        These are jobs that are not completed and not terminated in an error state
        """
        raise NotImplementedError

    def get_all(self):
        """
        return all jobs in queue in any state for this user
        """
        raise NotImplementedError

    def get_completed(self):
        """
        return list of completed jobs.  return jobs in same format as get_pending
        """
        return []

    def get_pending_count(self, pending_only=False):
        """
        number of pending jobs
        """
        queue = self.get_pending()
        if not pending_only:
            return len(queue)
        else:
            if self.PENDING_MARKER is None:
                return len(queue), 0  # no pending setup, act as if all are running
            else:
                return len(queue), len([j for j in queue if j[-1] == self.PENDING_MARKER])

    def get_errors(self):
        """
        jobs with errors, where sublists have ordering such as the first two items are [user,jobid,...]
        These are terminated jobs in an error state
        """
        raise NotImplementedError

    def get_error_count(self):
        """
        the number of errors
        """
        return len(self.get_errors())

    def get_subprocess_kwargs(self):
        return {'env': os.environ.copy(), 'shell': True}

class SlurmEngine(Engine):
    PENDING_MARKER = "PENDING"
    COMPLETED_MARKER = "COMPLETED"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.total_errors = 0
        self.last_check = 0
        self.last_queue = None

    def submit_job(self, path, name=None):
        current_dir = os.getcwd()
        working_dir = os.path.dirname(path)
        job_filename = os.path.basename(path)
        out = None
        try:
            os.chdir(working_dir)
            if name is None:
                name = job_filename
            cmd = "sbatch --job-name={name} {jobfile}".format(name=name, jobfile=job_filename)
            try:
                out = subprocess.check_output(
                    cmd, **self.get_subprocess_kwargs()).decode()
            except subprocess.CalledProcessError as e:
                self.total_errors += 1
                logging.warning("Problem (count={0}) calling sbatch: {1}".format(self.total_errors, e))
                if self.total_errors > MAX_ERRORS:
                    raise
                return None
        finally:
            os.chdir(current_dir)
        return out.split()[-1]

    def _check_queue(self, states, partition=None, qos=None, all_users=False):
        if time.time() > self.last_check + self.min_scan_interval:
            try:
                cmd = SACCT_CMD
                if not all_users:
                    cmd += ' -u %s' % getpass.getuser()
                queue = subprocess.check_output(
                    cmd, **self.get_subprocess_kwargs()).decode()
            except subprocess.CalledProcessError as e:
                try:
                    queue = subprocess.check_output(
                        SQUEUE_CMD, **self.get_subprocess_kwargs()).decode()
                except subprocess.CalledProcessError as e:
                    self.total_errors += 1
                    logging.warning("Problem (count={0}) calling squeue: {1}.  Using previous queue".format(self.total_errors, e))
                    if self.total_errors > MAX_ERRORS:
                        raise
            self.lines = [l.split() for l in str(queue).split("\n")]
            self.last_check = time.time()

        state_tuple = tuple(states.split(","))
        if all_users is False:
            username = getpass.getuser()
            self.last_queue = [line for line in self.lines[:-1] if line[0].startswith(username)]
        else:
            self.last_queue = self.lines[:-1]

        if partition is None and qos is None:
            return [job_line for job_line in self.last_queue if job_line[-1].startswith(state_tuple)]

        if partition is not None and qos is None:
            return [job_line for job_line in self.last_queue if job_line[-1].startswith(state_tuple) and partition in job_line[-3].split(',')]

        if partition is not None and qos is not None:
            return [job_line for job_line in self.last_queue if job_line[-1].startswith(state_tuple) and partition in job_line[-3].split(',') and qos == job_line[-2]]

        if partition is None and qos is not None:
            return [job_line for job_line in self.last_queue if job_line[-1].startswith(state_tuple) and qos == job_line[-2]]

    def update_job(self, jobid, partition, qos='normal', exc_node_list_string=None, nice=None):
        UPDATE_COMMAND = 'scontrol update job={0} partition={1} QOS={2}'.format(jobid, partition, qos)
        if nice is not None:
            UPDATE_COMMAND = UPDATE_COMMAND + ' nice={0}'.format(nice)
        if exc_node_list_string is not None:
            UPDATE_COMMAND = UPDATE_COMMAND + ' ExcNodeList={0}'.format(exc_node_list_string)
        try:
            subprocess.check_output(
                UPDATE_COMMAND, **self.get_subprocess_kwargs()).decode()
        except subprocess.CalledProcessError as e:
            self.total_errors += 1
            logging.warning("Problem (count={0}) calling scontrol update: {1}".format(self.total_errors, e))
            if self.total_errors > MAX_ERRORS:
                raise
            return None
        return None

    def get_completed(self):
        return self._check_queue(SLURM_COMPLETED_STATES)

    def get_pending(self):
        return self._check_queue(SLURM_RUN_STATES)

    def get_errors(self):
        return self._check_queue(SLURM_ERROR_STATES)

    def get_all(self):
        return self._check_queue(",".join([SLURM_COMPLETED_STATES, SLURM_RUN_STATES, SLURM_ERROR_STATES]))


class TorqueEngine(Engine):
    PENDING_MARKER = "Q,H"
    COMPLETED_MARKER = "C"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.total_errors = 0
        self.last_check = 0
        self.last_queue = None

    def submit_job(self, path, name=None):
        current_dir = os.getcwd()
        working_dir = os.path.dirname(path)
        job_filename = os.path.basename(path)
        out = None
        try:
            os.chdir(working_dir)
            if name is None:
                name = job_filename
            cmd = "qsub -N {name} {jobfile}".format(name=name, jobfile=job_filename)
            try:
                out = subprocess.check_output(
                    cmd, **self.get_subprocess_kwargs()).decode()
            except subprocess.CalledProcessError as e:
                self.total_errors += 1
                logging.warning("Problem (count={0}) calling sbatch: {1}".format(self.total_errors, e))
                if self.total_errors > MAX_ERRORS:
                    raise
                return None
        finally:
            os.chdir(current_dir)
        return out.split('.')[0]

    def _check_queue(self, states, queue=None, all_users=False):
        if time.time() > self.last_check + self.min_scan_interval:
            try:
                queued_calcs = subprocess.check_output(
                    QSTAT_CMD, **self.get_subprocess_kwargs()).decode()
            except subprocess.CalledProcessError as e:
                self.total_errors += 1
                logging.warning("Problem (count={0}) calling qstat: {1}.  Using previous queue".format(self.total_errors, e))
                if self.total_errors > MAX_ERRORS:
                    raise
            self.lines = [l.split() for l in str(queued_calcs).split("\n") if len(l) > 0 and l[0].isdigit()]
            self.last_check = time.time()

        state_tuple = tuple(states.split(","))
        if all_users is False:
            username = getpass.getuser()
            self.last_queue = [line for line in self.lines[:-1] if line[1].startswith(username)]
        else:
            self.last_queue = self.lines[:-1]

        if queue is None:
            return [job_line for job_line in self.last_queue if job_line[-2].startswith(state_tuple)]

        if queue is not None:
            return [job_line for job_line in self.last_queue if job_line[-2].startswith(state_tuple) and queue in job_line[2]]

    def get_completed(self):
        return self._check_queue(TORQUE_COMPLETED_STATES)

    def get_pending(self):
        return self._check_queue(TORQUE_RUN_STATES)

    def get_errors(self):
        return self._check_queue(TORQUE_ERROR_STATES)

    def get_all(self):
        return self._check_queue(",".join([TORQUE_COMPLETED_STATES, TORQUE_RUN_STATES, TORQUE_ERROR_STATES]))


class TestEngine(Engine):

    def submit_job(self, path, name=None):
        print("pretended to submit: " + path)

    def get_pending(self):
        return []

    def get_errors(self):
        return []


class LocalEngine(Engine):

    def __init__(self, *args, **kwargs):
        self.total_errors = 0
        self.errors = []
        self.current_job_id = 0

    def submit_job(self, path, name=None):
        current_dir = os.getcwd()
        working_dir = os.path.dirname(path)
        job_filename = os.path.basename(path)
        out = None
        self.current_job_id += 1
        try:
            os.chdir(working_dir)
            if name is None:
                name = job_filename
            cmd = "bash {jobfile}".format(jobfile=job_filename)
            try:
                out = subprocess.check_output(
                    cmd, stderr=subprocess.STDOUT,
                    **self.get_subprocess_kwargs()).decode()
                with open(os.path.join(working_dir, "job_output.txt"), 'w') as fp:
                    fp.write(out)
            except subprocess.CalledProcessError as e:
                self.total_errors += 1
                logging.warning("Error: {0} (count={1})\ncalling script: {2} ".format(e, self.total_errors, path))
                if self.total_errors > MAX_ERRORS:
                    raise
                return None
        finally:
            os.chdir(current_dir)
        return str(self.current_job_id)

    def get_pending(self):
        return []

    def get_all(self):
        return []

    def get_errors(self):
        username = getpass.getuser()
        return [(username, 0, e) for e in self.errors]
