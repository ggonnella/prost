#!/usr/bin/env bash

USAGE="Start the database server.\n"
USAGE+="\n"
USAGE+="Usage:\n  $0 <datadir> <dbname> <dbuser> <dbpass>\n"
USAGE+="\n"
USAGE+="Example:\n  $0 $HOME/prostdb_datadir prostdb prostuser prostpass\n"
USAGE+="\n"
USAGE+="Note:\n  The data directory must already been prepared\n"
USAGE+="  A database '<dbname>' must exist\n"
USAGE+="  An account '<dbuser>' with full privileges on '<dname>'\n"
USAGE+="  and password '<dbpass>' must exist.\n"
USAGE+="  This can be done e.g. using ./create_database.sh\n"

if [ $# -ne 4 ]; then
    echo -e $USAGE
    exit 1
fi

MSGPFX="[ $0 ] "

DATADIR=$1

if [ ! -d $DATADIR ]; then
    echo "${MSGPFX}Error: data directory $DATADIR does not exist."
    exit 1
fi

DBNAME=$2
DBUSER=$3
DBPASS=$4

INSTALL_ERRLOG=prost.db.server.log
INSTALL_SOCKET=prost.db.sock
INSTALL_PID=prost.db.server.pid

function start_server {
    echo -n "${MSGPFX}Step: start the database server... "
    TEMPFILE=$(mktemp)
    cmd="mysqld_safe --user=$USER \
                --datadir=$DATADIR \
                --pid-file=$DATADIR/$INSTALL_PID \
                --log-error=$DATADIR/$INSTALL_ERRLOG \
                --socket=$DATADIR/$INSTALL_SOCKET"
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
    while [ ! -e $DATADIR/$INSTALL_PID ]; do
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
echo "Values for the configuration file (config.yaml):"
echo ""
echo "  dbname: $DBNAME"
echo "  dbuser: $DBUSER"
echo "  dbpass: $DBPASS"
echo "  dbpath: $DATADIR"
echo "  dbsocket: $DATADIR/$INSTALL_SOCKET"
