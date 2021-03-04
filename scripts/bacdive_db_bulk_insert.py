#!/usr/bin/env python3
"""
Bulk insert Bacdive data.
The tables must exist already.

Usage:
  bacdive_db_bulk_insert.py [options] <dbuser> <dbpass> <dbname> <dbsocket> <file>

Arguments:
  dbuser:    database user to use
  dbpass:    password of the database user
  dbname:    database name
  dbsocket:  connection socket file
  file:      output of json_to_tabular.py

Options:
  --verbose, -v    be verbose
  --version, -V    show script version
  --help, -h       show this help message
"""
import MySQLdb
import tqdm
import sh
from docopt import docopt
from schema import Schema, And, Use, Or
import os
import json

def main(arguments):
  db = MySQLdb.connect(host="localhost",
                       user=arguments["<dbuser>"],
                       passwd=arguments["<dbpass>"],
                       db=arguments["<dbname>"],
                       unix_socket=arguments["<dbsocket>"],
                       use_unicode=True)
  cursor = db.cursor()
  noflines=int(str(sh.wc("-l", arguments["<file>"])).split(" ")[0])
  linedata = {}
  tablename = None
  with open(arguments["<file>"], 'r') as f:
    for line in tqdm.tqdm(f, total=noflines):
      line = line.rstrip()
      if line == "###":
        columns = linedata.keys()
        query = f"INSERT INTO {tablename}("
        query += ", ".join(columns)
        query += ") VALUES("
        query += ", ".join([f"%({k})s" for k in columns])
        query += ")"
        cursor.execute(query, linedata)
        linedata = {}
      else:
        elems = line.split("\t")
        tablename = elems[0]
        linedata[elems[1]] = json.loads(elems[2])
  cursor.close()
  db.commit()

def validated(arguments):
  schema = Schema({"<file>": open,
                   "<dbuser>": And(str, len),
                   "<dbpass>": And(str, len),
                   "<dbname>": And(str, len),
                   "<dbsocket>": And(str, len, os.path.exists),
                   str: object})
  return schema.validate(arguments)

if __name__ == "__main__":
  arguments = docopt(__doc__, version="0.1")
  main(validated(arguments))
