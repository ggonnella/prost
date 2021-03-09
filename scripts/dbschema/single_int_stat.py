#!/usr/bin/env python3
"""
Single int stat is a stat about a genome which is expressed by an integer
and for which there is a single number per genome.

Examples:
  - number of chromosomes
  - length of the genome
  - number of tRNA genes

A genome is thereby identified by its NCBI assembly accession
"""
from sqlalchemy.sql import func
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Sequence, Column, Integer, String, \
                       Index, DateTime, Date, Text, Boolean, Enum
from sqlalchemy.orm import relationship
from sqlalchemy_repr import PrettyRepresentableBase
import enum

Base = declarative_base(cls=PrettyRepresentableBase)

utf8_cs_args = {'mysql_charset': 'utf8', 'mysql_collate': 'utf8_bin'}

class SingleIntStat(Base):
  __tablename__ = 'single_int_stat'
  accession = Column(String(64), primary_key=True)
  statname = Column(String(128), primary_key=True)
  time_updated = Column(DateTime,
      server_default=func.now(), onupdate=func.now())
  value = Column(Integer)
  source = Column(String(512))
  __table_args__ = utf8_cs_args

