  # - environment
  #   - django intro
  #   - db credentials
  #   - git
  #     - forking
  # - updating
  #   - see docs/git.md

  # - molgen
  #  - structure of a
  # - scanconfigs
  # - automating jobs
  #   - requestjobs
  #   - buildjobs
  #   - parsejobs
  # - shell
  #   - filter
  #   - sort
  #   - csv dump
  # - web


# To start here, you have already run the environment setup

# make a place to create our new jobdirs
mkdir -p jobdirs


if [ "$(ls -A jobdirs)" ]; then
   echo "the directory, jobdirs, should be empty for this tutorial to run."
   exit
fi

./manage.py migrate

# first generate molecules - this makes 9!
./manage.py molgen example molgen/tests/fragments_small.json --recipe=test

# now load up the chemconfigs
./manage.py scanconfigs --update

# request some conformers on our new molecules
./manage.py requestjobs example conformer


# and now build them
./manage.py buildjobs example jobdirs

# Do the work!  (this is where you would use Odyssey or other muscle machine in real life)
# here, we just run the conformer generator.  This is a good example that you don't *always*
# need to run a job on the cluster
pushd jobdirs/*NSPMIYGKQJPBQR-TULZNQERNA-N
bash job.sh
popd

# now we parse those results back into the database
./manage.py parsejobs example jobdirs/

#Rinse, Repeat!
# now request a DFT to be run on the conformer (just one now)
# This time, we specify the parent config, this ensures that we don't request this job on *every* geom object in the DB
# right now, that's fine since we only have conformers in the db, but later, that's not what we want to do
./manage.py requestjobs example bp86_6-31gs_opt_qchem -p conformer

# just like before, build any unclaimed jobs
./manage.py buildjobs example jobdirs

# like any good cooking show, we've already run a job on odyssey ahead of time
cp tutorial/prebaked_dft/* jobdirs/*bp86_6-31gs_opt_qchem_NSPMIYGKQJPBQR-TULZNQERNA-N

# now parse the results again.
./manage.py parsejobs example jobdirs/