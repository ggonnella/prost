#!/usr/bin/env python3
"""
Compare a numerical value to a list of values

Usage:
  cmpstat_distribution.py [options] <value> <listfn>

Arguments:
  value:     value of interest
  listfn:    file with list to compare value of interest to, one value per line

Options:
  --verbose, -v    be verbose
  --version, -V    show script version
  --help, -h       show this help message
"""

from docopt import docopt
from schema import Schema, Or, Use, Optional
import scipy.stats

def values_from_listfile(f):
  line = f.readline().rstrip()
  try:
    result = [int(line)]
    useint = True
  except ValueError:
    result = [float(line)]
    useint = False
  if useint: result += [int(line.rstrip())   for line in f]
  else:      result += [float(line.rstrip()) for line in f]
  return result

def main(args):
  values = values_from_listfile(args["<listfn>"])
  value = args["<value>"]
  d = scipy.stats.describe(values)
  stdev = scipy.stats.tstd(values)
  zscoremin = (d.minmax[0]-d.mean)/stdev
  zscoremax = (d.minmax[1]-d.mean)/stdev
  zscore = (value-d.mean)/stdev
  percentile = scipy.stats.percentileofscore(values, value)
  print(f"Description of the group of values:")
  print(f"  Number of values in the group: {d.nobs}")
  print(f"  Range (min/avg/max): {d.minmax[0]}--{d.mean}--{d.minmax[1]}")
  print(f"  Standard deviation: {stdev}")
  print(f"  Z-score(min): {zscoremin}")
  print(f"  Z-score(max): {zscoremax}")
  print(f"Comparison to the value of interest:")
  print(f"  Z-score(value): {zscore}")
  print(f"  Z-score(value) - Z-score(min): {zscore-zscoremin}")
  print(f"  Z-score(max) - Z-score(value): {zscoremax-zscore}")
  print(f"  Percentile: {percentile}")
  interest = 0
  if value == d.mean:
    fact = "Value equal to the mean"
  elif value < d.mean:
    fact = "Value smaller than the mean"
    if value >= d.minmax[0]:
      if percentile < 25:
        interest = 1
        fact = "Value in the first quartile"
        if percentile < 10:
          interest = 2
          fact = "Value in the first decile"
    else:
      interest = 3
      fact = "Value smaller than any other"
      zdist = zscoremin-zscore
      if zdist > 1:
        interest = 4
        fact = "Value smaller than any other by more than 1 st.dev"
  else:
    fact = "Value larger than the mean"
    if value <= d.minmax[1]:
      if percentile > 75:
        interest = 1
        fact = "Value in the highest quartile"
        if percentile > 90:
          interest = 2
          fact = "Value in the highest decile"
    else:
      interest = 3
      fact = "Value larger than any other"
      zdist = zscore-zscoremax
      if zdist > 1:
        interest = 4
        fact = "Value larger than any other by more than 1 st.dev"
  print(f"Interest level: {interest}")
  print(f"Fact: {fact}")


def validated(args):
  schema = Schema({
    "<value>": Or(Use(int), Use(float)),
    "<listfn>": Use(open),
    Optional(str): object})
  return schema.validate(args)

if "snakemake" in globals():
  args = {
    "<value>": snakemake.params.value,
    "<listfn>": snakemake.input.list,
    }
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__, version="0.1")
  main(validated(args))
