#!/usr/bin/env python3
"""
Load data into DB from a TSV file.

Usage:
  db_load_data.py [options] <dbuser> <dbpass> <dbname> <dbsocket>
                            <tsv> <table> <columns> [<columns>...]

Arguments:
  dbuser:    database user to use
  dbpass:    password of the database user
  dbname:    database name
  dbsocket:  connection socket file
  tsv:       tsv file to load
  table:     table name
  columns:   either a list of columns or the path to a python module
             (see --dbschema)

Options:
  --ncbidmp, -n    use NCBI taxonomy dmp files delimiter settings, i.e.
                   <TAB>|<TAB> delimiter; <TAB>|<EOL> line ends
  --ignore, -i     use IGNORE on repeated primary key (default: REPLACE)
  --dbschema, -d   interpret the column argument as a python module name
                   which provides (1) a dictionary tablename2class, which
                   maps table names to classes and (2) a file_column_names()
                   class method for that class, which returns the list
                   of columns to use
  --set, -s        file with values to add to all records
                   (tsv with lines: column_name <TAB> value)
  --dropkeys, -k   drop non-unique indices before inserting
                   and re-compute them after inserting
  --verbose, -v    be verbose
  --version, -V    show script version
  --help, -h       show this help message
"""
import MySQLdb
from docopt import docopt
from schema import Schema, And, Or, Use, Optional
import os
import importlib

def column_names_from_dbschema(filename, tablename):
  spec = importlib.util.spec_from_file_location("dbschema", filename)
  dbschema = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(dbschema)
  klass = dbschema.tablename2class[tablename]
  return klass.file_column_names()

def load_data_set_from_file(fixed):
  result = ""
  if fixed:
    with open(fixed) as f:
      setelems = []
      for line in f:
        elems = line.rstrip().split("\t")
        setelems.append(f"{elems[0]} = \"{elems[1]}\" ")
      result += "SET "+", ".join(setelems)
  return result

def load_data_sql(datafile, table, columns, fixed, ignore,
                  dropkeys, ncbidmp):
  result = []
  result.append("SET foreign_key_checks = 0;")
  if dropkeys:
    result.append(f"ALTER TABLE {table} DISABLE KEYS;")
  sql = f"LOAD DATA LOCAL INFILE '{datafile}' "
  sql += "IGNORE " if ignore else "REPLACE "
  sql += f" INTO TABLE {table} "
  if ncbidmp:
    sql += r"FIELDS TERMINATED BY '\t|\t' "
    sql += r"LINES TERMINATED BY '\t|\n' "
  sql +="("+",".join(columns)+") "
  sql += load_data_set_from_file(fixed)
  sql += ";"
  result.append(sql)
  if dropkeys:
    result.append(f"ALTER TABLE {table} ENABLE KEYS;")
  result.append("SET foreign_key_checks = 1;")
  return result

def main(args):
  db = MySQLdb.connect(host="localhost",
                       user=args["<dbuser>"],
                       passwd=args["<dbpass>"],
                       db=args["<dbname>"],
                       unix_socket=args["<dbsocket>"],
                       use_unicode=True)
  cursor = db.cursor()
  if args["--dbschema"]:
    columns = column_names_from_dbschema(args["<columns>"][0],
                                         args["<table>"])
  else:
    columns = args["<columns>"]
  statements = load_data_sql(args["<tsv>"], args["<table>"], columns,
                      args["--set"], args["--ignore"],
                      args["--dropkeys"], args["--ncbidmp"])
  for statement in statements:
    cursor.execute(statement)
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
                   "--dbschema": Or(None, True, False),
                   "--ignore": Or(None, True, False),
                   "--dropkeys": Or(None, True, False),
                   "--ncbidmp": Or(None, True, False),
                   "--set": Or(None, open),
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
    "<columns>": snakemake.params.get("columns",
                   [snakemake.params.get("dbschema", None)]),
    "--dbschema": False,
    "--ignore": snakemake.params.get("ignore", False),
    "--dropkeys": snakemake.params.get("dropkeys", False),
    "--ncbidmp": snakemake.params.get("ncbidmp", False),
    "--set": snakemake.input.get("common_values", None)}
  if snakemake.params.get("dbschema", False):
    args["--dbschema"] = True
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__, version="0.1")
  main(validated(args))
