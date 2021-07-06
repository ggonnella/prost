#!/usr/bin/env python3
import math
import scipy.stats

values = []
with open("examples") as f:
  for line in f:
    values.append(
        [float(v) for v in line.split("#")[0].rstrip().split(" ")])
print(values)

LOG_BASE=2
MAX_SCORE=4

def entropy(freq, base=LOG_BASE):
  return scipy.stats.entropy(freq, base=base)

def max_entropy(freq, base=LOG_BASE):
  return math.log(len(freq), base)

def nonzero(freq):
  return [v for v in freq if v > 0]

def onezero(freq):
  return [v for v in freq if v > 0] + [0.0]

def max_entropy_nozero(freq, base=LOG_BASE):
  return max_entropy(nonzero(freq), base)

def ic(freq, base=LOG_BASE):
  return max_entropy(freq, base)-entropy(freq, base)

def norm_ic(freq, base=LOG_BASE):
  return 1-(entropy(freq, base)/max_entropy(freq, base))

def max_entropy_onezero(freq, base=LOG_BASE):
  return max_entropy(onezero(freq), base)

def complement(freq):
  return [1.0-x for x in freq]

def complementpow(freq, exp):
  return [pow(1.0-x, exp) for x in freq]

def scaledcompl(freq, base=LOG_BASE, maxscore=MAX_SCORE):
  return [x * maxscore for x in complementpow(freq, 4)]

from icecream import ic

for freq in values:
  print("----------------")
  ic(freq)
  ic(complement(freq))
  ic(complementpow(freq, 2))
  ic(complementpow(freq, 3))
  ic(entropy(freq))
  ic(ic(freq))
  ic(norm_ic(freq))
  ic(max_entropy(freq))
  ic(max_entropy_nozero(freq))
  ic(max_entropy_onezero(freq))
  ic(scaledcompl(freq))
