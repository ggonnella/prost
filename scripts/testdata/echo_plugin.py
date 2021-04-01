#!/usr/bin/env python3
"""
Echoes the input, for test purposes
"""

from pathlib import Path

ID =      "echo"
VERSION = "1.0"
INPUT = "anything"
OUTPUT =  ["echo"]

def compute(unit, **kwargs):
  return [unit], None

def compute_id(filename):
  return Path(filename).stem
