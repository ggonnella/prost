#!/usr/bin/env bash

USAGE="Start the database server.\n"
USAGE+="\n"
USAGE+="Usage:\n  $0 <datadir> <dbhost> <dbport> <dbsocket> <dbname> <dbuser> <dbpass>\n"
USAGE+="\n"
USAGE+="Example:\n"
USAGE+="  $0 $HOME/prostdb_datadir $HOME/prostdb_datadir/prost.db.socket\\\n"
USAGE+="                           localhost 3306 prostdb prostuser prostpass\n"
USAGE+="\n"
USAGE+="Note:\n  The data directory must already been prepared\n"
USAGE+="  A database '<dbname>' must exist\n"
USAGE+="  An account '<dbuser>' with full privileges on '<dname>'\n"
USAGE+="  and password '<dbpass>' must exist.\n"
USAGE+="  This can be done e.g. using ./create_database.sh\n"

if [ $# -ne 7 ]; then
    echo -e $USAGE
    exit 1
fi

MSGPFX="[ $0 ] "

DATADIR=$1
DBHOST=$2
DBPORT=$3
DBSOCKET=$4
DBNAME=$5
DBUSER=$6
DBPASS=$7

ERRLOG=$DBSOCKET.server.log
PID=$DBSOCKET.pid

function start_server {
    echo -n "${MSGPFX}Step: start the database server... "
    TEMPFILE=$(mktemp)
    cmd="mysqld_safe --user=$USER \
                --datadir=$DATADIR \
                --pid-file=$PID \
                --log-error=$ERRLOG \
                --socket=$DBSOCKET"
    ( ($cmd 2>&1) > $TEMPFILE ) &
    if [ $? -ne 0 ]; then
      echo -e "ERROR\n"
      echo -e "Command line: $cmd\n"
      echo -e "Command output:\n"
      cat $TEMPFILE
      rm $TEMPFILE
      exit 1
    fi
    MAX_TRIES=10
    while [ ! -e $PID ]; do
      sleep 1
      MAX_TRIES=$((MAX_TRIES-1))
      if [ $MAX_TRIES -eq 0 ]; then
        echo -e "ERROR\n"
        echo -e "Command line: $cmd\n"
        echo -e "Command output:\n"
        cat $TEMPFILE
        rm $TEMPFILE
        exit 1
      fi
    done
    rm $TEMPFILE
    echo "OK"
}

start_server
echo "Connection data for the configuration file:"
echo ""
echo "dbname: $DBNAME"
echo "dbuser: $DBUSER"
echo "dbpass: $DBPASS"
echo "dbsocket: $DBSOCKET"
echo "dbhost: $DBHOST"
echo "dbport: $DBPORT"
echo "dbpath: $DATADIR"
