#!/usr/bin/env python3
"""
DB Schema for tables representing per-genome statistics.
"""
from sqlalchemy.sql import func
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Sequence, Column, Integer, Float, String, DateTime
from sqlalchemy_repr import PrettyRepresentableBase
from sqlalchemy.dialects.mysql import DOUBLE

Base = declarative_base(cls=PrettyRepresentableBase)

utf8_cs_args = {'mysql_charset': 'utf8', 'mysql_collate': 'utf8_bin'}

class SingleIntStat(Base):
  """
  Single int stat is a stat about a genome which is expressed by an integer
  and for which there is a single number per genome.

  Examples:
    - number of chromosomes
    - length of the genome
    - number of tRNA genes

  A genome is thereby identified by its NCBI assembly accession
  """
  __tablename__ = 'single_int_stat'
  accession = Column(String(64), primary_key=True)
  statname = Column(String(128), primary_key=True)
  time_updated = Column(DateTime,
      server_default=func.now(), onupdate=func.now())
  value = Column(Integer)
  source = Column(String(512))
  __table_args__ = utf8_cs_args

class SingleFloatStat(Base):
  """
  Single float stat is a stat about a genome which is expressed by a float
  and for which there is a single number per genome.

  Examples:
    - GC percentage

  A genome is thereby identified by its NCBI assembly accession
  """
  __tablename__ = 'single_float_stat'
  accession = Column(String(64), primary_key=True)
  statname = Column(String(128), primary_key=True)
  time_updated = Column(DateTime,
      server_default=func.now(), onupdate=func.now())
  value = Column(DOUBLE)
  source = Column(String(512))
  __table_args__ = utf8_cs_args

tablename2class = {
  'single_int_stat': SingleIntStat,
  'single_float_stat': SingleFloatStat,
}

