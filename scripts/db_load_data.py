#!/usr/bin/env python3
"""
Load data into DB from a TSV file.

Usage:
  db_load_data.py [options] <dbuser> <dbpass> <dbname> <dbsocket> <tsv> <table> <columns> [<columns>...]

Arguments:
  dbuser:    database user to use
  dbpass:    password of the database user
  dbname:    database name
  dbsocket:  connection socket file
  tsv:       tsv file to load
  table:     table name
  columns:   list of columns

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
  query = "LOAD DATA LOCAL INFILE '{}' ".format(args["<tsv>"])
  query += "REPLACE INTO TABLE {} (".format(args["<table>"])
  query += ",".join(args["<columns>"])
  query += ");"
  cursor.execute(query)
  cursor.close()
  db.commit()

def validated(args):
  schema = Schema({"<dbuser>": And(str, len),
                   "<dbpass>": And(str, len),
                   "<dbname>": And(str, len),
                   "<dbsocket>": And(str, len, os.path.exists),
                   "<tsv>": And(str, open),
                   "<table>": And(str, len),
                   "<columns>": And(list, len),
                   Optional(str): object})
  return schema.validate(args)

if "snakemake" in globals():
  args = {
    "<dbuser>": snakemake.config["dbuser"],
    "<dbpass>": snakemake.config["dbpass"],
    "<dbname>": snakemake.config["dbname"],
    "<dbsocket>": snakemake.input.socket,
    "<tsv>": snakemake.input.tsv,
    "<table>": snakemake.params.table,
    "<columns>": snakemake.params.columns}
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__, version="0.1")
  main(validated(args))
