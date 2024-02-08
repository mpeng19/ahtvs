in order to make yohur own changes to this repo, you will need to fork it on bitbucket.

from then on, the process of updating will go as follows:

- one time

        git remote add upstream http://your_username@bitbucket.org/aspuru-guzik/a2g2.git

- every time to update

        git checkout master
        git fetch upstream
        git merge upstream/master

    - then add any new pip requirements (from environments/a2g2)
            
            pip install -r requirements.txt

    - and migrate the database for any updates to the model (from djangochem dir)
            
            ./manage.py migrate



