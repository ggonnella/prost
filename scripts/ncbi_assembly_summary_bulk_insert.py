#!/usr/bin/env python3
"""
Bulk insert NCBI assembly summary.
The table must exist already.

Usage:
  ncbi_taxonomy_db_bulk_insert.py [options] <dbuser> <dbpass> <dbname> <dbsocket> <file> <database> <domain>

Arguments:
  dbuser:    database user to use
  dbpass:    password of the database user
  dbname:    database name
  dbsocket:  connection socket file
  file:      filename
  database:  refseq or genbank
  domain:    bacteria or archaea

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
from schema import Schema, And, Use, Or
import os

def elems2values(arguments, elems):
  result = [elems[0],
            arguments["<database>"],
            arguments["<domain>"]]
  for i in range(len(elems)):
    if elems[i] == "na":
      elems[i] = None
    elif elems[i] == "":
      elems[i] = None
  for i in range(22-len(elems)):
    elems.append(None)
  result += elems[1:]
  return result

def main(arguments):
  db = MySQLdb.connect(host="localhost",
                       user=arguments["<dbuser>"],
                       passwd=arguments["<dbpass>"],
                       db=arguments["<dbname>"],
                       unix_socket=arguments["<dbsocket>"],
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
  if arguments["--update"]:
    query +=" ON DUPLICATE KEY UPDATE "
    query += ", ".join(f"{k}=VALUES({k})" for k in columns[1:])
  batch = []
  noflines=int(str(sh.wc("-l", arguments["<file>"])).split(" ")[0])
  with open(arguments["<file>"], 'rb') as f:
    i = 0
    for line in tqdm.tqdm(f, total=noflines):
      i += 1
      line = line.decode('cp1252', errors="ignore")
      if line[0] != "#":
        if i % arguments["--batch"] == 0:
          cursor.executemany(query, batch)
          batch.clear()
        elems = line.rstrip().split("\t")
        batch.append(elems2values(arguments, elems))
  if len(batch):
    cursor.executemany(query, batch)
  cursor.close()
  db.commit()

def validated(arguments):
  schema = Schema({"<file>": open,
                   "<dbuser>": And(str, len),
                   "<dbpass>": And(str, len),
                   "<dbname>": And(str, len),
                   "<dbsocket>": And(str, len, os.path.exists),
                   "<database>": Or("refseq", "genbank"),
                   "<domain>": Or("bacteria", "archaea"),
                   "--batch": And(Use(int), lambda n: n>0),
                   str: object})
  return schema.validate(arguments)

if __name__ == "__main__":
  arguments = docopt(__doc__, version="0.1")
  main(validated(arguments))
