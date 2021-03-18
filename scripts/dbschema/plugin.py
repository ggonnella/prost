#!/usr/bin/env python3
"""
Plugins for the computation of assembly attributes.
"""
from sqlalchemy.sql import func
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Sequence, Column, Integer, Float, String, DateTime,\
                       Boolean, Text, Enum
from sqlalchemy_repr import PrettyRepresentableBase
from sqlalchemy.dialects.mysql import DOUBLE

Base = declarative_base(cls=PrettyRepresentableBase)

utf8_cs_args = {'mysql_charset': 'utf8', 'mysql_collate': 'utf8_bin'}

class Plugin(Base):
  """
  Describes a plugin for assembly attributes computation.
  """
  __tablename__ = "plugin"
  id = Column(String(256), primary_key=True)
  version = Column(String(64), primary_key=True)
  input = Column(String(512), nullable=False)
  output = Column(String(512), nullable=False)
  method = Column(Text(4096))
  implementation = Column((Text(4096)))
  parameters = Column((Text(4096)))
  req_software = Column(Text(4096))
  req_hardware = Column(Text(4096))
  advice = Column(Text(4096))
  __table_args__ = utf8_cs_args

"""
- id:              fas_stats_nim
  version:         0.1.0
  input:           genome sequence; Fasta; optionally Gzip-compressed
  output:          genome_size, GC_content
  method:          count bases
  implementation:  NIM module, lookup table
  parameters:      gzipped (Boolean)
  req_software:    nim >= 1.2.8; nimble:zip >= 0.3.1
  req_hardware:    NULL
  advice:          preferred to fasgz_stats_posix, if nim is available,
                   since it is faster

- id:              fas_stats_posix
  version:         0.1.0
  input:           genome sequence; Fasta; optionally Gzip-compressed
  output:          genome_size, GC_content
  method:          count bases
  implementation:  posix tools pipe using sh library
  comparison:      NULL
  req_software:    grep, tr, wc, zcat
  req_hardware:    NULL
  parameters:      gzipped (Boolean)

- id:              refseq_gff_count
  version:         0.1.0
  input:           genome annotation; from Refseq; gff, optionally gzipped
  output:          feature_type_count
  method:          - counts features of a predefined set of feature_type
                   - feature_type is defined as in GFF (i.e. term from the SOFA
                     subset of the SO ontology)
                   - to avoid counting portions of the same feature multiple
                     times, children of the same type of the same gene
                     are counted only once
                   - for some feature_type (e.g. rRNA), a more specific
                     feature_type (e.g. rRNA_16S) is assigned using information
                     from the attributes (e.g. product_name)
  implementation:  NULL
  comparison:      NULL
  req_software:    NULL
  req_hardware:    NULL
  parameters:      NULL
"""
