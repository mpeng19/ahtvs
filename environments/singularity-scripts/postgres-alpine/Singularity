BootStrap: docker
From: postgres:alpine

%startscript
/docker-entrypoint.sh postgres -h $HOSTNAME

%environment
export HOSTNAME=localhost
if [ -f /postgresrc ]; then 
    . /postgresrc 
fi

