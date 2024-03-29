#!/usr/bin/env python3
"""
Bulk insert Bacdive data.
The tables must exist already.

Usage:
  bacdive_db_bulk_insert.py [options] {db_args_usage} <data>

Arguments:
{db_args}
  file:      output of json_to_tabular.py

Options:
{common}
"""
import MySQLdb
import tqdm
import sh
from lib import db, scripts
import json
import snacli

def main(args):
  db = MySQLdb.connect(host="localhost",
                       user=args["<dbuser>"],
                       passwd=args["<dbpass>"],
                       db=args["<dbname>"],
                       unix_socket=args["<dbsocket>"],
                       use_unicode=True)
  cursor = db.cursor()
  noflines=int(str(sh.wc("-l", args["<data>"])).split(" ")[0])
  linedata = {}
  tablename = None
  with open(args["<data>"], 'r') as f:
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

def validated(args):
  return scripts.validate(args, db.args_schema, {"<data>": open})

with snacli.args(db.snake_args,
                 docvars={"common": scripts.args_doc,
                          "db_args": db.args_doc,
                          "db_args_usage": db.args_usage},
                 input=["<data>"], version="0.1") as args:
  if args: main(validated(args))
