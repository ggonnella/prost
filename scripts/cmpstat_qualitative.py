#!/usr/bin/env python3
"""
Compare the value of a categorical variable to a list of values

Usage:
  cmpstat_categories.py [options] <value> <listfn>

Arguments:
  value:     value of interest
  listfn:    file with list to compare value of interest to, one value per line

Options:
  --bin BINS   distribute numeric values into bins:
               given the argument N0,N1,...,Ni,...,Nn, the i-th bin is:
                 - i == 0:     [value < N0]
                 - 0 < i < n:  [N(i-1) <= value < Ni]
                 - i == n:     [value  >= Nn]
  --kmeans K   use k-means clustering to assign numeric values to K bins
{common}
"""

from docopt import docopt
from schema import Or, Use, And
from lib import snake, scripts
from collections import Counter
from bisect import bisect_right

def get_values(value, f):
  values = [value]
  values += [line.rstrip() for line in f]
  try:
    int_values = [int(v) for v in values]
    return int_values[0], int_values[1:]
  except ValueError:
    try:
      float_values = [float(v) for v in values]
      return float_values[0], float_values[1:]
    except ValueError:
      return value, values[1:]

def compare(value, values):
  interest = 0
  print(f"Number of elements of the group: {len(values)}")
  if value not in values:
    print("Value not present in the group")
    if len(values) > 10:
      interest = 4
    elif len(values) > 4:
      interest = 3
    elif len(values) > 2:
      interest = 2
    else:
      interest = 1
  else:
    counts = Counter(values)
    freq_before = {k: c/(len(values)) for k, c in counts.items()}
    counts[value] += 1
    freq_after = {k: c/(len(values)+1) for k, c in counts.items()}
    print("Value present in the group")
    print(f"Frequency of the value in the group: {freq_before[value]*100:.1f}%")
    print(f"Frequency of the value in (group & value): {freq_after[value]*100:.1f}%")
    if freq_after[value] < 0.1:
      interest = 3
    elif freq_after[value] < 0.25:
      interest = 2
    elif freq_after[value] < 0.5:
      interest = 1
  print(f"Interest level: {interest}")

def bin_values(value, values, bin_thresholds):
  if isinstance(value, str):
    raise ValueError("Cannot bin non-numeric values")
  return bisect_right(bin_thresholds, value), \
         [bisect_right(bin_thresholds, v) for v in values]

def kmeans_fit_values(value, values, k):
  from sklearn.cluster import KMeans
  import numpy as np
  a = np.array([value] + values).reshape(-1, 1)
  kmeans = KMeans(n_clusters=k).fit_predict(a)
  return kmeans[0], list(kmeans[1:])

def main(args):
  value, values = get_values(args["<value>"], args["<listfn>"])
  if args["--bin"]:
    value, values = bin_values(value, values, args["--bin"])
  elif args["--kmeans"]:
    value, values = kmeans_fit_values(value, values, args["--kmeans"])
  compare(value, values)

def parse_thresholds(s):
  return sorted([float(t) for t in set(s.split(","))])

def validated(args):
  return scripts.validate(args, {"<value>": str,
                                 "<listfn>": Use(open),
                                 "--bin": Or(None, Use(parse_thresholds)),
                                 "--kmeans": Or(None, And(Use(int),
                                                  lambda x: x > 1))})

if "snakemake" in globals():
  args = snake.args(snakemake, params = ["<value>", "--bin"],
                    input=[("<listfn>", "list")])
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__.format(common=scripts.args_doc), version="0.1")
  main(validated(args))
