#!/usr/bin/env bash

USAGE="Create a database and an user for interacting with that database.\n"
USAGE+="\n"
USAGE+="Usage:\n  $0 <datadir> <password> <dbname> <dbuser> <dbpass>\n"
USAGE+="\n"
USAGE+="Example:\n  $0 $HOME/prostdb_datadir prostdb prostuser prostpass\n"
USAGE+="\n"
USAGE+="Note:\n  The data directory must already been prepared\n"
USAGE+="  A full privileges account '$USER' must exist, with password '<password>'.\n"
USAGE+="  This can be done e.g. using ./install_database.sh\n"

if [ $# -ne 5 ]; then
    echo -e $USAGE
    exit 1
fi

MSGPFX="[ $0 ] "

DATADIR=$1

if [ ! -d $DATADIR ]; then
    echo "${MSGPFX}Error: data directory $DATADIR does not exist."
    exit 1
fi

PASSWORD=$2
DBNAME=$3
DBUSER=$4
DBPASS=$5

INSTALL_ERRLOG=create_database.server.log
INSTALL_SOCKET=create_database.server.socket
INSTALL_PID=create_database.server.pid

function check_command { cmdname=$1; shift; helpmsg=$*;
  echo -n "${MSGPFX}Checking that $cmdname is available... "
  which $cmdname &> /dev/null
  if [ $? -ne 0 ]; then
    echo "ERROR"
    echo "Error: $cmdname not found."
    echo $helpmsg
    exit 1
  fi
  echo "OK"
}

function run_command { step=$1; shift; cmd=$*;
    echo -n "${MSGPFX}Step: $step... "
    TEMPFILE=$(mktemp)
    ($cmd 2>&1) > $TEMPFILE
    if [ $? -ne 0 ]; then
      echo -e "ERROR\n"
      echo -e "Command line: $cmd\n"
      echo -e "Command output:\n"
      cat $TEMPFILE
      rm $TEMPFILE
      exit 1
    fi
    rm $TEMPFILE
    echo "OK"
}

function run_sql { step=$1; shift; cmd=$*;
    echo -n "${MSGPFX}Step: $step... "
    TEMPFILE=$(mktemp)
    (mysql -u$USER --password=$PASSWORD \
            -S$DATADIR/$INSTALL_SOCKET \
            --verbose --execute="$cmd;" 2>&1) > $TEMPFILE
    if [ $? -ne 0 ]; then
      echo -e "ERROR\n"
      echo -e "SQL: $cmd\n"
      echo -e "Command output:\n"
      cat $TEMPFILE
      rm $TEMPFILE
      exit 1
    fi
    rm $TEMPFILE
    echo "OK"
}

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
run_sql "create database '$DBNAME'" \
  "CREATE DATABASE \`$DBNAME\`"
#run_sql "create user '$DBUSER'" \
#  "CREATE USER '$DBUSER' IDENTIFIED BY '$DBPASS'"
run_sql "grant usage to '$DBUSER'" \
  "GRANT USAGE ON *.* TO '$DBUSER'@localhost IDENTIFIED BY '$DBPASS'; "
run_sql "grant privileges to '$DBUSER'" \
  "GRANT ALL PRIVILEGES ON \`$DBNAME\`.* TO '$DBUSER'@localhost; FLUSH PRIVILEGES"

run_command "shutting down the database server" \
  "mysqladmin -u $USER --socket=$DATADIR/$INSTALL_SOCKET \
              --password=$PASSWORD shutdown"
