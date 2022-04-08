#!/usr/bin/env python3
"""
Download the pages with IDs lists from Bacdive
using the REST API and extract from it the list
of IDs.

As the download is very slow, take an argument,
which is the first page from which to start downloading
and process only pages after it.

Usage:
  bacdive_download_id_list.py [options] <username> <password> <idlist>

Arguments:
  username        username on bacdive
  password        password on bacdive
  idlist          file containing the IDs
                  (new IDs are appended if file exists)

Options:
  --done VALUE    last ID list page already processed
{common}
"""
import requests
from requests.auth import HTTPBasicAuth
from schema import Use, Or
from lib import scripts
import sys
import re
import snacli

def get(url, headers, credentials, verbose):
  if verbose:
    sys.stderr.write(f"# GET: {url}\n")
  response = requests.get(url, headers = headers, auth=credentials)
  if verbose:
    sys.stderr.write(f"# response status code: {response.status_code}\n")
  if response.status_code != 200:
    response.raise_for_status()
  return response.json()

def compute_first_url(baseurl, args, headers, credentials):
  prevurl = baseurl
  url = None
  if args["--done"] is not None:
    # determine in this case the "next" to be processed
    # by reading it from the last processed one
    donepage = args["--done"]
    if donepage > 0:
      prevurl = f"{prevurl}?page={donepage}"
    json = get(prevurl, headers, credentials, args["--verbose"])
    url = json['next']
    if args["--verbose"]:
      sys.stderr.write("# Accessing last processed URL to determine next URL\n")
    if url is None:
      if args["--verbose"]:
        sys.stderr.write("# No next URL to process\n")
        sys.stderr.write("# Everything already processed\n")
    else:
      if args["--verbose"]:
        sys.stderr.write("# Next URL is: {}\n".format(url))
  return prevurl, url

def main(args):
  if args["--verbose"]:
    sys.stderr.write("# This script updates the list of IDs from Bacdive\n")
    sys.stderr.write(f"# Username: {args['<username>']}\n")
    sys.stderr.write(f"# Password: {args['<password>']}\n")

  headers = {'Accept': 'application/json'}
  baseurl='https://bacdive.dsmz.de/api/bacdive/bacdive_id/'
  credentials = HTTPBasicAuth(args["<username>"], args["<password>"])
  prevurl, url = compute_first_url(baseurl, args, headers, credentials)

  if url:
    while url:
      json = get(url, headers, credentials, args["--verbose"])
      if args["--verbose"]:
        sys.stderr.write("# Processing results...\n")
      for result in json['results']:
        if args["--verbose"]:
          sys.stderr.write("## Result: \n"+ result['url'])
        args["<idlist>"].write(result['url'].split('/')[-2]+"\n")
      prevurl = url
      url = json['next']
      if args["--verbose"]:
        if url is None:
          sys.stderr.write("# No next URL to process\n")
          sys.stderr.write("# Everything processed\n")
        else:
          sys.stderr.write("# Next URL is: {}\n".format(url))
  m = re.search(r"\d+", prevurl)
  print(m.group(0))

def validated(args):
  return scripts.validate(args, {
    "<idlist>": Use(lambda fn: open(fn, "a+")),
    "<username>": str, "<password>": str,
    "--done": Or(None, int)})

with snacli.args(log=["<idlist>"],
                 params=["<username>", "<password>", "--done"],
                 docvars={"common": scripts.args_doc},
                 version="0.1") as args:
  if args: main(validated(args))
