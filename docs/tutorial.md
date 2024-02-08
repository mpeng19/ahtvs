Outline

  - environment
    - django intro
    - db credentials
    - git
      - forking
  - updating
    - see docs/git.md
  - molgen
  - scanconfigs
  - automating jobs
    - requestjobs
    - buildjobs
    - parsejobs
  - shell
    - filter
    - sort
    - csv dump
  - web

10. Setting up your environment
    - http://bitbucket.com/aspuru-guzik/a2g2
    - click source -> environments -> a2g2
    - follow instructions in README
      - istall Miniconda (python3)
      - git clone a2g2 repo
      - cd a2g2/environments/
      - run setup_env.sh
      - quick conda intro (see more [here](http://conda.pydata.org/docs/intro.html))
      - most important commands:
        - `source activate a2g2` to enter the environment (note you will see the name of the env in parens at your prompt)
        - `source deactivate` to exit the environment (or just close your terminal)
      - install chemaxon (this installs the tools you need to render images of your structures), see: docs/chemaxon.md
      - install homebrew and obabel

        brew install open-babel

15. setup database


    - ./manage.py migrate
    - ./manage.py createsuperuser  - this isn't the database user, this is the admin user for the website that you may or may not actually use.

   - django settings file intro
     - how to set it
     - environment variables
     - default.py
     - how manage.py refers to it with --settings=

16. Molgen
    - ./manage.py molgen example molgen/tests/fragments_small.json --recipe=test

28. Run your first job
    - (one time) ./manage.py scanconfigs --update
    - ./manage.py requestjobs example conformer
    - ./manage.py buildjobs example jobdirs

30. Executing jobdirs
    - getting work to odyssey
      - sshfs
      - scp
    - job_manager.py
    - the batch trick

50. djangochem
    - django project
    - why a web tool, django?
      	  - free settings managements - picking databases, db backup,
	  - free command line tool setup
	  - ./manage.py help

60. chemconfigs
    - link to jinja2


65. REST api
    - ./manage.py createsuperuser
    - ./manage runserver
    - and visit http://127.0.0.1:8000
    - in the api, click the dropdown menu

70. Appendix / seperate topics
    - sharing code and forking repos
    - overview of a2g2 directory

