======
djangochem
======

djangochem is a django app for the storage of theoretical chemistry calculations and data


setting up test environment:
docker run -p 27017:27017 --name some-mongo -d mongo
docker run --name some-postgres -p 5432:5432 -e POSTGRES_PASSWORD=postgres -d postgres


requirements:

docker (it's still in beta for mac but it should be out soon) - 
    This is really just for running the databases.  You can also setup the databases in other ways, or use dbs that are hosted on other machines.  Docker is just super easy.
    
python3


run tests
---------

./manage.py test
