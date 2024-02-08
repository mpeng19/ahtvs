Environment:

There are 2 ways to get going.  You should choose the docker way unless you have reason not to.

----------------
## the EASY Docker way!
1. Install docker. https://docs.docker.com/engine/installation/

   As of July 28, 2016, docker for mac and windows is out of beta.  That means, don't install virtualbox, boot2docker, or docker-machine.  All that is gone.  If you find docs online referencing those tools, walk away.
   
   On Linux, you would need to install `docker-compose`. Make
   sure your compose version is > 1.11.
   https://docs.docker.com/compose/install/

2. Clone this repo
HTTPS (you need to create token/password on GitLab.com): `git clone https://gitlab.com/dennisaag/a2g2.git`
SSH: `git clone git@gitlab.com:dennisaag/a2g2.git`

3. Run docker and install environments.
  * `cd environments/a2g2-docker-dev/`
  * `./run.sh`.  
  * This will spin up a postgres database and install your whole environment and drop you at a the command line inside the docker where you run all the django commands that comprise the a2g2 tool set.  

4. If the database were not setup before, run `django-admin migrate` to setup the database. 

5. Testing. Run `django-admin test`. If this is your first run, 

6. From that docker command line run `bash tutorial/auto_tutorial.sh`  This will instantly run a tutorial.  You are going to want to come back and run this line by line.  if you want to re-run line by line, use `./manage.py flush` to erase the whole local database and start again.  and `rm -r jobdirs` to erase all the job directories created by the tutorial.

7. Assuming you don't yet reset the database as shown in step 5, you can now run `./run_notebook_docker.sh`. This will fire up a jupyter notebook.  Now point your browser to http://localhost:9999 (with token) and you are about ready to do some science! 

**Not updated**
## The HARD WAY.  If you want this environment directly installed in your dev box.  
1. install miniconda from here: http://conda.pydata.org/miniconda.html
    it doesn't matter if you pick the python2 or 3 release because both can create environments of the other type.  If in doubt, use the python3.  This code base is going that way.

1.5.  Clone this repo!

2. import the provided conda environment and enter it.  In this case, I am calling the environemnt 'a2g2', but feel free to use your own naming convention.  This will be a python 3 environment.  run these commands from this directory.  You can just paste all the lines.  Assumed current working dir is the same as this README


        ./setup_env.sh

3. From now on you can enter the environment with `source activate a2g2` and leave the environment with `source deactivate`.

4. Install Database.

    - download and install postgres http://postgresapp.com/

    - (alternate) If you prefer homebrew, run `brew install postgres`
    - brew tap homebrew/services
    - https://github.com/Homebrew/homebrew-services

    createuser postgres -P
    (provide postgres as password)

    - Linux Postgres Install:
        https://www.digitalocean.com/community/tutorials/how-to-install-and-use-postgresql-on-ubuntu-14-04

5. Install a GUI database browser:  https://eggerapps.at/postico/

6. Please also read docs/chemaxon.md and follow those instructions to install chemaxon if you want a nice way to see your molecules.  It is also required for hosting a local website to browse your molecule database.

7. Optional

        brew tap mcs07/cheminformatics
        brew install open-babel --with-python
        echo 'import site; site.addsitedir("/usr/local/lib/python2.7/site-packages")' >> ~/miniconda3/envs/a2g2/lib/python3.5/site-packages/homebrew.pth
    - see more here: https://github.com/mcs07/homebrew-cheminformatics
    
    
    
##Docker
###Docker for development (local a2g2 code)

1. Pull the a2g2 code
    HTTPS: `git clone https://gitlab.com/dennisaag/a2g2.git`  (you need to make a token on gitlab to use https)
    SSH: `git clone git@gitlab.com:dennisaag/a2g2.git`
1. 
    `cd a2g2/environments/a2g2-docker-dev/`

1. Build the images first time from scratch

    `sudo docker-compose build --pull` or `sudo docker-compose build --pull --no-cache`  
    
1. Run the docker and enter the shell:

    `sudo docker-compose run a2g2`

1. Test the code:
    `django-admin test`

1. To remove dockers
    1. without removing database file: `sudo docker-compose down`

    1. remove docker and volume database
    `sudo docker-compose down -v`


### Docker with a2g2 codebase
TBD




 
