#!/usr/bin/env python3
"""
Load data into DB from a TSV file.

Usage:
  db_load_data.py [options] {db_args_usage}
                            <tsv> <table> <columns> [<columns>...]

Arguments:
{db_args}
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
{common}
"""
from docopt import docopt
from schema import And, Or, Use
from lib import scripts, snake, db, mysql, sqlwriter

def main(args):
  headerpfx = args["--headerpfx"] if args["--skipheader"] else ""
  columns = args["<columns>"][0] if args["--dbschema"] else args["<columns>"]
  statements = sqlwriter.load_data_sql(args["<tsv>"], args["<table>"], columns,
                      args["--skipfields"], args["--set"], args["--ignore"],
                      args["--dropkeys"], args["--ncbidmp"], headerpfx)
  mysql.connect_and_execute(args, statements)

def validated(args):
  return scripts.validate(args, db.args_schema,
                  {"<tsv>": And(str, open),
                   "<table>": And(str, len),
                   "<columns>": And(list, len),
                   "--dbschema": Or(None, True, False),
                   "--ignore": Or(None, True, False),
                   "--dropkeys": Or(None, True, False),
                   "--skipfields": Or(And(None, Use(lambda n: [])),
                     And(str, Use(lambda l: [int(e) for e in l.split(",")]),
                         lambda l: all(e > 0 for e in l))),
                   "--ncbidmp": Or(None, True, False),
                   "--skipheader": Or(None, True, False),
                   "--headerpfx": Or(And(None, Use("#")), And(str, len)),
                   "--set": Or(None, open)})

if "snakemake" in globals():
  args = snake.args(snakemake, db.snake_args,
        input=["<tsv>", ("--set", "common_values")],
        params=["<table>", "<columns>", "--dbschema", "--ignore",
                "--dropkeys", "--ncbidmp", "--skipheader", "--headerpfx",
                "--skipfields"])
  if args["--dbschema"]:
    args["<columns>"] = args["--dbschema"]
    args["--dbschema"] = True
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__.format(db_args = db.args_doc, common = scripts.args_doc,
           db_args_usage = db.args_usage), version="0.1")
  main(validated(args))
