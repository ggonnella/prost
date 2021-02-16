#!/usr/bin/env python3
"""
Download the pages with IDs lists from Bacdive
using the REST API and extract from it the list
of IDs.

As the download is very slow, take an argument,
which is the first page from which to start downloading
and process only pages after it.

Usage:
  download_id_list.py [options] <username> <password> <idlist>

Arguments:
  username        username on bacdive
  password        password on bacdive
  idlist          file containing the IDs
                  (new IDs are appended if file exists)

Options:
  --done VALUE    last ID list page already processed
  --verbose, -v   be verbose
  --help, -h      show this help message
  --version, -V   show the script version
"""
import requests
from requests.auth import HTTPBasicAuth
from docopt import docopt
import sys
import re

def get(url, headers, credentials, verbose):
  if verbose:
    sys.stderr.write(f"# GET: {url}\n")
  response = requests.get(url, headers = headers, auth=credentials)
  if verbose:
    sys.stderr.write(f"# response status code: {response.status_code}\n")
  if response.status_code != 200:
    response.raise_for_status()
  return response.json()

def compute_first_url(baseurl, arguments, headers, credentials):
  prevurl = baseurl
  url = None
  if arguments["--done"]:
    # determine in this case the "next" to be processed
    # by reading it from the last processed one
    donepage = int(arguments["--done"])
    if donepage > 0:
      prevurl = f"{prevurl}?page={donepage}"
    json = get(prevurl, headers, credentials, arguments["--verbose"])
    url = json['next']
    if arguments["--verbose"]:
      sys.stderr.write("# Accessing last processed URL to determine next URL\n")
    if url is None:
      if arguments["--verbose"]:
        sys.stderr.write("# No next URL to process\n")
        sys.stderr.write("# Everything already processed\n")
    else:
      if arguments["--verbose"]:
        sys.stderr.write("# Next URL is: {}\n".format(url))
  return prevurl, url

def main(arguments):
  if arguments["--verbose"]:
    sys.stderr.write("# This script updates the list of IDs from Bacdive\n")
    sys.stderr.write("# IDs will be stored in {}\n".format(
                        arguments["<idlist>"]))
    sys.stderr.write(f"# Username: {arguments['<username>']}\n")
    sys.stderr.write(f"# Password: {arguments['<password>']}\n")

  headers = {'Accept': 'application/json'}
  baseurl='https://bacdive.dsmz.de/api/bacdive/bacdive_id/'
  credentials = HTTPBasicAuth(arguments["<username>"], arguments["<password>"])
  prevurl, url = compute_first_url(baseurl, arguments, headers, credentials)

  if url:
    with open(arguments["<idlist>"], "a+") as f:
      while url:
        json = get(url, headers, credentials, arguments["--verbose"])
        if arguments["--verbose"]:
          sys.stderr.write("# Processing results...\n")
        for result in json['results']:
          if arguments["--verbose"]:
            sys.stderr.write("## Result: \n"+ result['url'])
          f.write(result['url'].split('/')[-2]+"\n")
        prevurl = url
        url = json['next']
        if arguments["--verbose"]:
          if url is None:
            sys.stderr.write("# No next URL to process\n")
            sys.stderr.write("# Everything processed\n")
          else:
            sys.stderr.write("# Next URL is: {}\n".format(url))
  m = re.search(r"\d+", prevurl)
  print(m.group(0))

if __name__ == "__main__":
  arguments = docopt(__doc__, version="0.1")
  arguments["--verbose"]=True
  main(arguments)
