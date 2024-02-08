a2g2 environment setup
======================

Setting up a2g2 environment by `shub://verysure/rdkit-django` and `shub://verysure/postgres-alpine`.
Include easy-to-use `a2g2` wrapper script.

Development
-----------

The environments depend on singularity containers: `verysure/rdkit-django` and `verysure/postgres-alpine`.
Clone the github repository to build the singularity containers from the recipes. A backup copy of the contents
are also included in `a2g2/environments/singularity-scripts/`. Pull request to the changes are welcomed.
If you would like to maintain the singularity build repository, send an email to "ttttonywu [at] gmail [dot] com" 
to gain administrative access of the github repository.

Installation
------------
1. Install `singularity`, see [link](https://www.sylabs.io/guides/3.0/user-guide/quick_start.html#quick-installation-steps). Unfortunately, newest version only supports linux.
2. Go to `a2g2/environments/a2g2-singularity`.
3. Run `bash install.sh` or `./install.sh` install the environments. (optional) Edit `ENV_PATH` in `install.sh` to change default installation path.
4. Edit `/path/to/install/settings`. Please change your `DB_NAME` and `POSTGRES_PASSWORD`. (your installed path, not the one in a2g2 repository).
5. (optional) Set `a2g2` alias in `~/.bashrc` for quick access.

Usage
-----
1. Start postgres: `/path/to/install/a2g2 pg start` to start database service. (must)
2. Run `/path/to/install/a2g2 test` to test your a2g2 setup.
3. Run `/path/to/install/a2g2 shell` to execute commands in a2g2 environment.
4. Run `/path/to/install/a2g2 jupyter` (in your choice of notebook paths) to start jupyter server.
5. Run `/path/to/install/a2g2` to see all available commands.

Database
--------
1. You can change and modify the database with `a2g2`. 
2. You can create a new database by `/path/to/install/a2g2 db create <db_name>`.

