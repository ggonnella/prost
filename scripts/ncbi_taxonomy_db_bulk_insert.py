#!/usr/bin/env python3
"""
Bulk insert NCBI taxonomy data from dump files.
The table must exist already.

Usage:
  ncbi_taxonomy_db_bulk_insert.py [options] <dbuser> <dbpass> <dbname> <dbsocket> <file> <table>

Arguments:
  dbuser:    database user to use
  dbpass:    password of the database user
  dbname:    database name
  dbsocket:  connection socket file
  file:      filename
  table:     table name

Options:
  --update, -u     update records if not unique (default: fail if not unique)
  --batch=N        number of records to pass to executemany [default: 20000]
  --verbose, -v    be verbose
  --version, -V    show script version
  --help, -h       show this help message
"""
import MySQLdb
import tqdm
import sh
from docopt import docopt
from schema import Schema, And, Use, Optional
from dbschema.ncbi_taxonomy_db import tablename2class
import os

def main(args):
  db = MySQLdb.connect(host="localhost",
                       user=args["<dbuser>"],
                       passwd=args["<dbpass>"],
                       db=args["<dbname>"],
                       unix_socket=args["<dbsocket>"],
                       use_unicode=True)
  columns = tablename2class[args["<table>"]].file_column_names()
  cursor = db.cursor()
  query = "INSERT INTO "
  query += args["<table>"] + "("
  query += ", ".join(columns)
  query += ") VALUES("
  query += ", ".join(["%s"] * len(columns))
  query += ")"
  if args["--update"]:
    query +=" ON DUPLICATE KEY UPDATE "
    query += ", ".join(f"{k}=VALUES({k})" for k in columns)
  print(query)
  batch = []
  noflines=int(str(sh.wc("-l", args["<file>"])).split(" ")[0])
  #if not args["--update"]:
  #  cursor.execute(f"ALTER TABLE {args['<table>']} DISABLE KEYS;")
  with open(args["<file>"]) as f:
    i = 0
    for line in tqdm.tqdm(f, total=noflines):
      i += 1
      if i % args["--batch"] == 0:
        cursor.executemany(query, batch)
        batch.clear()
      elems = line.rstrip()[:-2].split("\t|\t")
      batch.append(elems)
  if len(batch):
    cursor.executemany(query, batch)
  #if not args["--update"]:
  #  cursor.execute(f"ALTER TABLE {args['<table>']} ENABLE KEYS;")
  cursor.close()
  db.commit()

def validated(args):
  schema = Schema({"<file>": open,
                   "<dbuser>": And(str, len),
                   "<dbpass>": And(str, len),
                   "<dbname>": And(str, len),
                   "<dbsocket>": And(str, len, os.path.exists),
                   "<table>": And(str, len),
                   "--batch": And(Use(int), lambda n: n>0),
                   Optional(str): object})
  return schema.validate(args)

if "snakemake" in globals():
  args = {
    "<file>": snakemake.input.dump,
    "<dbuser>": snakemake.config["dbuser"],
    "<dbpass>": snakemake.config["dbpass"],
    "<dbname>": snakemake.config["dbname"],
    "<dbsocket>": snakemake.input.socket,
    "<table>": snakemake.params.table,
    "--batch": 20000}
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__, version="0.1")
  main(validated(args))
