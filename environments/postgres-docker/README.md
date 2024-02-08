Setting up Environment
======================

Set environment variables:

`export POSTGRES_PASSWORD=secret_adasda`

`export DB_NAME=a2g2_test`

`export DJANGO_SETTINGS_MODULE=djangochem.settings.postgres_docker`

(you can also set these in conda environment to autoload when sourcing a2g2_test, 
https://conda.io/docs/user-guide/tasks/manage-environments.html#saving-environment-variables)

Run postgres docker in background
`sudo -E docker-compose up -d`

Shutdown postgres docker
`sudo -E docker-compose down`

Shutdown postgres docker with removing database
`sudo -E docker-compose down -v`

Make new db:
`sudo docker exec -it a2g2_postgres createdb -U postgres a2g2_test`

`sudo docker exec -it a2g2_postgres createdb -U postgres $DB_NAME`

First time use: `django-admin migrate`


Import db data:

`psql -d a2g2_test -h localhost -p 5440 -U postgres -f /path/dump.sql`

Export db data:


For use in terminal:

export POSTGRES_PASSWORD=secret_adasda
export DB_NAME=a2g2_test
export DJANGO_SETTINGS_MODULE=djangochem.settings.postgres_docker
sudo -E docker-compose up -d
sudo docker exec -it a2g2_postgres createdb -U postgres a2g2_test

