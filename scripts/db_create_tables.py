#!/usr/bin/env python3
"""
Create all tables from a SqlAlchemy schema file.

Usage:
  db_create_tables.py [options] {db_args_usage} <schema> [<table>]

Arguments:
{db_args}
  schema:    python file containing the SqlAlchemy classes
  table:     create only specified table (optional, default: create all)

Requirements:
  the file must define a global variable Base from declarative_base() e.g.
  Base = declarative_base(cls=PrettyRepresentableBase)

Options:
  --drop, -d       drop table if exists (requires <table> argument)
{common}
"""
from sqlalchemy.orm.session import sessionmaker
from sqlalchemy import create_engine
from docopt import docopt
from schema import Or
from lib import db, snake, scripts
import importlib

def main(args):
  engine = create_engine(db.connstr_from(args), echo=args["--verbose"])
  spec = importlib.util.spec_from_file_location("models", args["<schema>"])
  models = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(models)
  if args["<table>"]:
    if args["--drop"]:
      query = "DROP TABLE IF EXISTS `{}`;".format(args["<table>"])
      Session = sessionmaker(bind=engine)
      session = Session()
      session.execute(query)
    models.Base.metadata.tables[args["<table>"]].create(bind=engine.engine)
  else:
    models.Base.metadata.create_all(engine.engine)

def validated(args):
  return scripts.validate(args, db.args_schema, {"<schema>": open,
    "<table>" : Or(None, str), "--drop": Or(None, True, False)})

if "snakemake" in globals():
  args = snake.args(snakemake, db.snake_args, input=["<schema>"],
                    params=["<table>", "--drop"])
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__.format(db_args = db.args_doc,
    db_args_usage = db.args_usage, common = scripts.args_doc), version="0.1")
  main(validated(args))
