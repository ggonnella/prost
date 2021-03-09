#!/usr/bin/env python3
"""
Bulk insert NCBI taxonomy data from dump files.
The table must exist already.

Usage:
  db_drop_table.py [options] <dbuser> <dbpass> <dbname> <dbsocket> <table>

Arguments:
  dbuser:    database user to use
  dbpass:    password of the database user
  dbname:    database name
  dbsocket:  connection socket file
  table:     table name

Options:
  --verbose, -v    be verbose
  --version, -V    show script version
  --help, -h       show this help message
"""
import MySQLdb
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
  cursor = db.cursor()
  query = "DROP TABLE `{}`;".format(args["<table>"])
  cursor.execute(query)
  cursor.close()
  db.commit()

def validated(args):
  schema = Schema({"<dbuser>": And(str, len),
                   "<dbpass>": And(str, len),
                   "<dbname>": And(str, len),
                   "<dbsocket>": And(str, len, os.path.exists),
                   "<table>": And(str, len),
                   Optional(str): object})
  return schema.validate(args)

if "snakemake" in globals():
  args = {
    "<dbuser>": snakemake.config["dbuser"],
    "<dbpass>": snakemake.config["dbpass"],
    "<dbname>": snakemake.config["dbname"],
    "<dbsocket>": snakemake.input.socket,
    "<table>": snakemake.params.table}
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__, version="0.1")
  main(validated(args))
