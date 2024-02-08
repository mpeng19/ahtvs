#Workers


A worker:

10. scans the database for a work_order that matches its own name
20. creates a job directory from that parent database document using template files
30. submits that directory to the job_manager
40. scans completed jobs folders
50. creates calculation documents from the resulting files using a loader module
60. saves those new documents in the database


###Configs

Each worker is specified by a json config file.  Config files live with their dependent files in subdirectories of the configs directory.  There are also database config files in the configs directory.  Inside a config subdirectory is everything to create jobs and return data to the database.

 - **templates** (for the jinja2 library) are used to create any files needed by the calculation software to run.
 - **loaders** are python modules that open files in the resulting directory and parse them to create database documents.

 Expect to see these files in a config dir:

10. `job.sh` - a template file to submit to slurm
20. `loader.py` - a python module that contains a function load_calc_list that accepts the path to the completed job directory, and the parent database document
30. `__init__.py` - these dirs must be python packages so we can import loader
40. `config.json` - a json formatted master config file that points to all these other template and loader files and tells the worker paths to use for storing job dirs: inbox, completed, archive, and error.
50. (optional) - extra template files such as a `.xyz` template and an `.inp` template.