
# job dict special keys
WORKER_NAME_KEY = "worker_name"
PARENT_KEY_KEY = "parent_key"
PARENT_CLASS_KEY = "parent_class"
JOB_CLASS_KEY = "job_class"
JOB_KEY = "job_key"
PROJECT_KEY = "project_name"
BATCH_KEY = "batch"

class LockError(Exception):
    pass

class JobNotFoundError(Exception):
    pass

class DatabaseInterface(object):
    def get_new_jobs(self, worker_name, project_name, min_priority=0, max_priority=10, limit=1):
        """
        returns calculation query objects containing calcs that this worker can work on
        """
        raise NotImplementedError

    def claim_job(self, jobspec):
        """ claiming a job prevents other interfaces from allowing it get claimed"""
        raise NotImplementedError

    def clear_job(self, jobspec):
        """ clear claimed a job"""
        raise NotImplementedError

    def mark_job_complete(self, jobinfo, worker_name):
        raise NotImplementedError

    def release_job(self, jobinfo, errors=False, ignore_lock=False):
        raise NotImplementedError

    def add_result(self, result, jobinfo, suffix=None):
        raise NotImplementedError
