#!/usr/bin/env python3
"""
Download the data for strains from Bacdive
using the REST API, output them as TSV.

Columns in the output TSV:
  (1) Bacdive ID
  (2) MD5 Hash of the downloaded JSON data
  (3) Timestamp of the download time
  (4) JSON data downloaded

Usage:
  bacdive_download_strain_info.py [options] <username> <password> <idlist>

Arguments:
  username       username on bacdive
  password       password on bacdive
  idlist         file containing the IDs to download

Options:
  --previous, -p FILENAME  text file containing previously downloaded data
{common}
"""

import requests
from requests.auth import HTTPBasicAuth
import sys
import hashlib
from lib import scripts
from datetime import datetime
from schema import Use, Or
import snacli

def main(args):
  if args["--verbose"]:
    sys.stderr.write("# This script downloads strain data from Bacdive\n")
    sys.stderr.write(f"# Username: {args['<username>']}\n")
    sys.stderr.write(f"# Password: {args['<password>']}\n")
  processed_ids = set()
  if args["--previous"]:
    for line in args["--previous"]:
      bacdiveid = line.split("\t")[0]
      processed_ids.add(bacdiveid)
  if args["--verbose"]:
    sys.stderr.write(\
        f"# N. previously downloaded strains: {len(processed_ids)}\n")
  headers = {'Accept': 'application/json'}
  credentials = HTTPBasicAuth(args["<username>"], args["<password>"])
  baseurl='https://bacdive.dsmz.de/api/bacdive/bacdive_id/'
  n_new_ids = 0
  n_old_ids = 0
  for line in args["<idlist>"]:
    bacdive_id = line.rstrip()
    if bacdive_id in processed_ids:
      n_old_ids += 1
    else:
      n_new_ids += 1
      url = baseurl + bacdive_id
      if args["--verbose"]:
        sys.stderr.write("# GET: "+url+"\n")
      response = requests.get(url, headers = headers, auth=credentials)
      if response.status_code != 200:
        response.raise_for_status()
      md5 = hashlib.md5(response.text.encode('utf-8')).hexdigest()
      ts = str(datetime.timestamp(datetime.now()))
      print("\t".join([bacdive_id, md5, ts, response.text]))
  if args["--verbose"]:
    sys.stderr.write("# All IDs processed\n")
    sys.stderr.write(f"# N. IDs ignored as already downloaded: {n_old_ids}\n")
    sys.stderr.write(f"# N. IDs for which info was downloaded: {n_new_ids}\n")

def validated(args):
  return scripts.validate(args, {"<idlist>": Use(open), "<username>": str,
    "<password>": str, "--previous": Or(None, Use(open))})

with snacli.args(input=["<idlist>"],
                 params=["<username>", "<password>", "--previous"],
                 docvars={"common": scripts.args_doc},
                 version="0.1") as args:
  if args: main(validated(args))
