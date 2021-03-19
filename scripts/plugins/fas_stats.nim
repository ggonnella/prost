##
## Compute simple sequence statistics for a Fasta file
##

import strutils, os
import math
import nimpy
import zip/gzipfiles
import tables

const bufsize = 2 ^ 16

type fas_state = enum indesc, newline, inseq

const isalphatab = block:
  var tmp: array[256, uint8]
  for i in 0 .. 255:
    tmp[i] = if isalphaascii(chr(i)): 1 else: 0
  tmp

const isgctab = block:
  var tmp: array[256, uint8]
  for i in 0 .. 255:
    tmp[i] = if chr(i) in ['G', 'C', 'g', 'c']: 1 else: 0
  tmp

template process_buffer(buffer: array[bufsize, char], nread: int,
                    count: var uint, state: var fas_state,
                    tab: array[256, uint8]) =
  for i in 0 ..< nread:
    var char = buffer[i]
    if state == indesc:
      if char == '\n':
        inc(state)
    elif state == newline:
      if char == '>':
        dec(state)
      else:
        inc(state)
        count += tab[ord(char)]
    else:
      if char == '\n':
        dec(state)
      else:
        count += tab[ord(char)]

template process_buffer2(buffer: array[bufsize, char], nread: int,
                    count_a: var uint,
                    count_b: var uint,
                    state: var fas_state,
                    tab_a: array[256, uint8],
                    tab_b: array[256, uint8]) =
  for i in 0 ..< nread:
    var char = buffer[i]
    if state == indesc:
      if char == '\n':
        inc(state)
    elif state == newline:
      if char == '>':
        dec(state)
      else:
        inc(state)
        count_a += tab_a[ord(char)]
        count_b += tab_b[ord(char)]
    else:
      if char == '\n':
        dec(state)
      else:
        count_a += tab_a[ord(char)]
        count_b += tab_b[ord(char)]

proc process_gzipped(fasfile: string, tab: array[256, uint8]): uint =
  var
    f = newGzFileStream(fasfile)
    buffer: array[bufsize, char]
    state: fas_state = indesc
  while true:
    var nread = f.readData(addr buffer, bufsize)
    if nread == 0:
      break
    process_buffer(buffer, nread, result, state, tab)
  f.close()

proc process_gzipped2(fasfile: string, tab_a: array[256, uint8],
                      tab_b: array[256, uint8]): tuple[a:uint, b:uint] =
  var
    f = newGzFileStream(fasfile)
    buffer: array[bufsize, char]
    state: fas_state = indesc
  while true:
    var nread = f.readData(addr buffer, bufsize)
    if nread == 0:
      break
    process_buffer2(buffer, nread, result.a, result.b, state, tab_a, tab_b)
  f.close()

proc process_file(fasfile: string, tab: array[256, uint8]): uint =
  var
    f: File
    buffer: array[bufsize, char]
  if open(f, fasfile):
    var state: fas_state = indesc
    while true:
      var nread = f.read_buffer(addr buffer, bufsize)
      if nread == 0:
        break
      process_buffer(buffer, nread, result, state, tab)
    f.close()

proc process_file2(fasfile: string, tab_a: array[256, uint8],
                   tab_b: array[256, uint8]): tuple[a:uint, b:uint] =
  var
    f: File
    buffer: array[bufsize, char]
  if open(f, fasfile):
    var state: fas_state = indesc
    while true:
      var nread = f.read_buffer(addr buffer, bufsize)
      if nread == 0:
        break
      process_buffer2(buffer, nread, result.a, result.b, state, tab_a, tab_b)
    f.close()

proc totallen(fasfile: string, gzipped: bool): uint {.exportpy.} =
  if gzipped: process_gzipped(fasfile, isalphatab)
  else: process_file(fasfile, isalphatab)

proc stats(fasfile: string, gzipped: bool): (int, float) {.exportpy.} =
  let (gc, total) = block:
    if gzipped: process_gzipped2(fasfile, isgctab, isalphatab)
    else: process_file2(fasfile, isgctab, isalphatab)
  return (total.int, gc.int/total.int)

proc analyze(filename: string,
             kwargs: Table[string, PyObject]):
              tuple[results: Table[string, Table[string, string]],
                    logs: Table[string, seq[string]]] {.exportpy.} =
  ## id: fas_stats_nim
  ## version: 0.1.0
  ## input: genome sequence; Fasta; optionally Gzip-compressed
  ## output: sequence_attribute.(genome_size,GC_content)
  ## method: count bases
  ## implementation: NIM module, lookup table
  ## parameters: force_uncompressed (bool)
  ## req_software: nim >= 1.2.8; nimble:zip >= 0.3.1
  ## advice: >
  ##   preferred to fasgz_stats_posix, if nim is available, since it is faster
  let force_uncompressed = block:
    if kwargs.has_key("force_uncompressed"):
      kwargs["force_uncompressed"].to(bool)
    else: false
  let (total, gc) = stats(filename, not force_uncompressed)
  result.results = {"sequence_attribute":
              {"genome_size": $(total),
               "gc_content": $(gc),
               "force": $(force_uncompressed),
               }.to_table}.to_table

proc gcpercent(fasfile: string, gzipped: bool): float {.exportpy.} =
  return stats(fasfile, gzipped)[1]

