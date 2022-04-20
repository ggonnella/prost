#!/usr/bin/env bash

USAGE="Prepare a database data directory and setup the database administration password.\n"
USAGE+="\n"
USAGE+="Usage:\n  $0 <data dir> <database password>\n"
USAGE+="\n"
USAGE+="Example:\n  $0 $HOME/prostdb_datadir mypassword\n"
USAGE+="\n"
USAGE+="Note:\n  The data directory must not exist yet.\n"

if [ $# -ne 2 ]; then
    echo -e $USAGE
    exit 1
fi

MSGPFX="[ $0 ] "

DATADIR=$1

if [ -d $DATADIR ]; then
    echo "${MSGPFX}Error: data directory $DATADIR already exists."
    exit 1
fi

PASSWORD=$2

if [ -z "$PASSWORD" ]; then
    echo "${MSGPFX}Error: password is empty."
    exit 1
fi

INSTALL_ERRLOG=install_database.server.log
INSTALL_SOCKET=install_database.server.socket
INSTALL_PID=install_database.server.pid

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

check_command mysql_install_db "Please install MariaDB first."
run_command "prepare the data directory" \
      "mysql_install_db --datadir=$DATADIR --user=$USER"

start_server
run_command "setting the password for database user '$USER'" \
  "mysqladmin -u $USER --socket=$DATADIR/$INSTALL_SOCKET password $PASSWORD"
run_command "shutting down the database server" \
  "mysqladmin -u $USER --socket=$DATADIR/$INSTALL_SOCKET \
              --password=$PASSWORD shutdown"

