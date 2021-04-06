#!/usr/bin/env python3
"""
Bulk insert NCBI assembly summary.
The table must exist already.

Usage:
  ncbi_taxonomy_db_bulk_insert.py [options] {db_args_usage}
                                  <file> <database> <domain>

Arguments:
{db_args}
  file:      filename
  database:  refseq or genbank
  domain:    bacteria or archaea

Options:
  --update, -u     update records if not unique (default: fail if not unique)
  --batch=N        number of records to pass to executemany [default: 20000]
{common}
"""
import MySQLdb
import tqdm
import sh
from docopt import docopt
from schema import And, Use, Or
from lib import db, snake, scripts

def elems2values(args, elems):
  result = [elems[0],
            args["<database>"],
            args["<domain>"]]
  for i in range(len(elems)):
    if elems[i] == "na":
      elems[i] = None
    elif elems[i] == "":
      elems[i] = None
  for i in range(22-len(elems)):
    elems.append(None)
  result += elems[1:]
  return result

def main(args):
  db = MySQLdb.connect(host="localhost",
                       user=args["<dbuser>"],
                       passwd=args["<dbpass>"],
                       db=args["<dbname>"],
                       unix_socket=args["<dbsocket>"],
                       use_unicode=True)
  columns = ["accession", "seqdb", "domain",
             "bioproject", "biosample", "wgs_master",
             "refseq_category", "taxid", "species_taxid",
             "organism_name", "infraspecific_name",
             "isolate", "version_status", "assembly_level",
             "release_type", "genome_rep", "seq_rel_date",
             "asm_name", "submitter", "gbrs_paired_asm",
             "paired_asm_comp", "ftp_path", "excluded_from_refseq",
             "relation_to_type_material"]
  cursor = db.cursor()
  query = "INSERT INTO ncbi_assembly_summary ("
  query += ", ".join(columns)
  query += ") VALUES("
  query += ", ".join(["%s"] * len(columns))
  query += ")"
  if args["--update"]:
    query +=" ON DUPLICATE KEY UPDATE "
    query += ", ".join(f"{k}=VALUES({k})" for k in columns[1:])
  batch = []
  noflines=int(str(sh.wc("-l", args["<file>"])).split(" ")[0])
  with open(args["<file>"], 'rb') as f:
    i = 0
    for line in tqdm.tqdm(f, total=noflines):
      i += 1
      line = line.decode('cp1252', errors="ignore")
      if line[0] != "#":
        if i % args["--batch"] == 0:
          cursor.executemany(query, batch)
          batch.clear()
        elems = line.rstrip().split("\t")
        batch.append(elems2values(args, elems))
  if len(batch):
    cursor.executemany(query, batch)
  cursor.close()
  db.commit()

def validated(args):
  return scripts.validate(args, db.args_schema,
          {"<file>": open, "<database>": Or("refseq", "genbank"),
           "<domain>": Or("bacteria", "archaea"),
           "--batch": Or(And(None, Use(lambda n: 20000)),
                                 And(Use(int), lambda n: n>0))})

if "snakemake" in globals():
  args = snake.args(snakemake, db.snake_args,
      input = [("<file>", "datasrc")],
      wildcards = [("<database>", "db"), "domain"],
      params = ["--update", "--batch"])
  main(args)
elif __name__ == "__main__":
  args = docopt(__doc__.format(db_args = db.args_doc,
     db_args_usage = db.args_usage, common = scripts.args_doc), version="0.1")
  main(validated(args))
