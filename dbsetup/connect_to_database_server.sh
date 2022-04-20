#!/usr/bin/env bash

USAGE="Interactive prompt to the database server.\n"
USAGE+="\n"
USAGE+="Usage:\n  $0 <datadir> <dbname> <dbuser> <dbpass>\n"
USAGE+="\n"
USAGE+="Example:\n  $0 $HOME/prostdb_datadir prostdb prostuser prostpass\n"
USAGE+="\n"
USAGE+="Note:\n  The database server must already be running\n"
USAGE+="  This can be done e.g. using ./start_database_server.sh\n"

if [ $# -ne 4 ]; then
    echo -e $USAGE
    exit 1
fi

MSGPFX="[ $0 ] "

DATADIR=$1
DBSOCK=$DATADIR/prost.db.sock

if [ ! -e $DBSOCK ]; then
    echo "${MSGPFX}Error: connection socket $DBSOCK does not exist."
    exit 1
fi

DBNAME=$2
DBUSER=$3
DBPASS=$4

mysql --socket=$DBSOCK --database=$DBNAME -u $DBUSER --password=$DBPASS

