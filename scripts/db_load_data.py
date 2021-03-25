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
  --skipheader     skip any line at the beginning of the file starting
                   with the header prefix (see --headerpfx);
                   internal lines starting with the prefix are not skipped
                   (consider preprocessing with grep -v instead in this case)
  --headerpfx PFX  header line prefix used by --skipheader (default: #)
  --skipfields L   comma-sep list of 1-based field numbers to skip
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
from docopt import docopt
from schema import Schema, And, Or, Use, Optional
import os
from lib import snake, mysql, sqlwriter

def main(args):
  headerpfx = args["--headerpfx"] if args["--skipheader"] else ""
  columns = args["<columns>"][0] if args["--dbschema"] else args["<columns>"]
  statements = sqlwriter.load_data_sql(args["<tsv>"], args["<table>"], columns,
                      args["--skipfields"], args["--set"], args["--ignore"],
                      args["--dropkeys"], args["--ncbidmp"], headerpfx)
  mysql.connect_and_execute(args, statements)

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
                   "--skipfields": Or(And(None, Use(lambda n: [])),
                     And(str, Use(lambda l: [int(e) for e in l.split(",")]),
                         lambda l: all(e > 0 for e in l)))
                   "--ncbidmp": Or(None, True, False),
                   "--skipheader": Or(None, True, False),
                   "--headerpfx": Or(And(None, Use("#")), And(str, len)),
                   "--set": Or(None, open),
                   Optional(str): object})
  return schema.validate(args)

if "snakemake" in globals():
  args = snake.args(snakemake,
        config=["<dbuser>", "<dbpass>", "<dbname>"],
        input=["<dbsocket>", "<tsv>", ("--set", "common_values")],
        params=["<table>", "<columns>", "--dbschema", "--ignore",
                "--dropkeys", "--ncbidmp", "--skipheader", "--headerpfx",
                "--skipfields"])
  if args["--dbschema"]:
    args["<columns>"] = args["--dbschema"]
    args["--dbschema"] = True
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__, version="0.1")
  main(validated(args))
