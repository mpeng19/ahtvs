Singularity Container
=====================

Singularity containter for postgres-alpine docker image ([link](https://hub.docker.com/_/postgres/)).
This singularity recipe modifies from the image to work in service mode. 
It also includes an optional `/postgresrc` to pass system environment variables.


Usage
-----

* To start the service, `singularity instance.start -B host_folder:/var/lib/postgresql/data postgres.img pg-database`.
* To pass variables for postgres database, use bind a file to `/postgresrc`.  
  For example, `-B yourrc:/postgresrc`. Some common variables include `PGPORT` and `HOSTNAME`.
  (Due to some reason, the default postgres starts in `0.0.0.0`. For security reason, this recipe
  uses `HOSTNAME` which defaults to `localhost`.)
* Use `SINGULARITY_BINDPATH='host_folder:container_folder,host_file:container_file'` for easier binding.



Building Manually
-----------------

To build the image, run `sudo singularity build <name.img> Singularity`.
See [Singularity](https://singularity.lbl.gov/) 
for more info. 


Acknowledge
-----------

The recipes build from many open source projects, including
* [Singularity](https://singularity.lbl.gov/)
* [postgres](https://www.postgresql.org/)
* [Alpine Linux](https://www.alpinelinux.org/)

