# TODO Items

July 2016 - This document outlines some good potential tasks for upcoming a2g2 maintainers from the perspective of the original author (tim).


- ### merging tool

    There a variety of situations where someone would want to merge the data from one database into another, potentially changing the associated project name (technically, changing the associated Group model).  For example, at the end of a project or a someone's postdoc, all their calculations should ideally be merged into a master group database.  This DB can then become a great resource for students who can benefit from a bunch of calculations, eg, machine learning or statistics on methods.

    - see https://docs.djangoproject.com/en/1.9/topics/db/multi-db/ for info on managing multiple databases.
    - It seems like this tool could be another command in the dbmover app (which current;y has a tool for migrating data from an old mongo db into the new postgres db)


- ### javascript datatables tool for the REST backend

    I have a sloppily written but functional example of this in a file called table.html in web/templates/  (found in either a2g2 or samsung repos.  same file.)  Two important doc sources:
    The current REST api is incomplete and lives in a22g/djangochem/restapi
    The main missing part is the ability to filter and sort.  It is browsable now if you `./manage.py runserver`

        - https://datatables.net/
        - http://www.django-rest-framework.org/

    You may see in the table.html example that I leverage the format=svg capability of the rest API to render nice images inside the tables.  the svg stuff is working in the new api.