#!/usr/bin/env python3
"""
Bulk insert NCBI taxonomy data from dump files.
The table must exist already.

Usage:
  db_drop_table.py [options] {db_args_usage} <table>

Arguments:
{db_args}
  table:     table name

Options:
{common}
"""
import MySQLdb
from docopt import docopt
from schema import And
from lib import snake, scripts, db

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
  return scripts.validate(args, db.args_schema, {"<table>": And(str, len)})

if "snakemake" in globals():
  args = snake.args(snakemake, db.snake_args, params = ["<table>"])
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__.format(db_args = db.args_doc, common = scripts.args_doc,
    db_args_usage = db.args_usage), version="0.1")
  main(validated(args))
