## Ferrocene Pipeline Example

This example demonstrates how to generate jobs for doing analysis on ferrocenes.

The `<my_*>` notation indicates that you should replace this with your settings/files.

### 2. Generate CP Fragments and Ferrocene Molecules via `molgen`

1. `./manage.py molgen -r ferrocene <my_project> <my_cp_fragments.json>`

### 3. Generate conformers of CP Fragments via `conformer`

1. Request conformer jobs on fragments:
  `./manage.py requestjobs <my_project> conformer -t ferrocene_fragment`

2. Build job dirs:
  `./manage.py buildjobs <my_project> <my_jobs_inbox>`

3. Run jobs (may take a bit)
  `cwd=$(pwd); for jobdir in $(find <my_jobs_inbox> -maxdepth 1); do cd $jobdir; bash job.sh; cd $cwd; done;`

4. Parse jobs
  `./manage.py parsejobs <my_project> <my_jobs_inbox> -r <my_jobs_archive>`

### 4. Generate Ferrocene conformers.
1. Request ferrocene_conformer job on ferrocene_molecule.
  `./manage.py requestjobs <my_project> ferrocene_conformer -t ferrocene`

2. Build job dirs:
  `./manage.py buildjobs <my_project> <my_jobs_inbox>`

3. Run jobs (may take a bit)
  `cwd=$(pwd); for jobdir in $(find <my_jobs_inbox> -maxdepth 1); do cd $jobdir; bash job.sh; cd $cwd; done;`

4. Parse jobs
  `./manage.py parsejobs <my_project> <my_jobs_inbox> -r <my_jobs_archive>`

At this point ferrocene conformers will exist in the database as `Geom` objects, with parentjob `ferrocene_conformer`.

Note that there can be *a lot* of conformers. It may make sense to do a first-pass dedupe/energy filtering before generating DFT jobs.

### 5. Build jobs for DFTs
  `./manage.py requestjobs <my_project> bp86_6-31gs_opt_qchem -p ferrocene_conformer`
The '`-p ferrocene_conformer`' bit is what tells the `requestjobs` command to only operate on the ferrocene conformers that we generated above.

At this point you can run jobs on odyssey as you normally would (probably via job_manager).
