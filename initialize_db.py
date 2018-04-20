#import os,sys
from sqlalchemy import Column,ForeignKey,Integer,String,Text,Enum,Float
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy import create_engine

Base = declarative_base()

class Organism(Base):
    __tablename__ = 'organisms'
    name = Column(Text,nullable=False)
    id = Column(Integer,primary_key=True,nullable=False)
    cds = relationship("CDS",back_populates='seqid')
    hmmHits = relationship("HMMhits",back_populates='hmmhit')
    rnaHits = relationship("RNAHits",back_populates='rnahit')


class CDS(Base):
    __tablename__= 'Seqs'
    seqid = Column(Integer,primary_key=True,nullable=False)
    orgname = Column(Text,ForeignKey('organisms.name'))
    gene = Column(Text)
    description = Column(Text)
    source = Column(Text)
    loc_start = Column(Integer)
    loc_end = Column(Integer)
    loc_strand = Column(Enum('+','-'))
    acc = Column(Text)
    lastscan = Column(Integer)
    naseq = Column(Text)
    aaseq = Column(Text,nullable=False)

    hmmHits = relationship("HMMhits",back_populates='hmmhit')

class HmmHits(Base):
    __tablename__='HMMhits'
    hmmhit = Column(Text,primary_key=True,nullable=False)
    orgname = Column(Text,ForeignKey('organisms.name'))
    seqid = Column(Integer,ForeignKey('Seqs.seqid'))
    hmmstart = Column(Integer)
    hmmend = Column(Integer)
    hmmlen = Column(Integer)
    geneStart = Column(Integer)
    geneEnd = Column(Integer)
    genelen = Column(Integer)
    evalue = Column(Float)
    score = Column(Float)
    bias = Column(Float)
    iscore = Column(Float)
    flags = Column(Integer)
    hmmcov = Column(Float)
    genecov = Column(Float)
    hmmLib = Column(Text)
    cds = relationship("Seqs",back_populates='orgname')

class RnaHits(Base):
    __tablename__ = 'RNAHits'
    rnahit = Column(Text)
    orgname = Column(Text, ForeignKey('organisms.name'))
    seqid = Column(Integer,primary_key=True)
    score = Column(Float)
    loc_start = Column(Integer)
    loc_end = Column(Integer)
    loc_strand = Column(Enum('+','-'))
    seq = Column(Text, nullable=False)

engine = create_engine('sqlite:///mlst_sql.db')
Base.metadata.create_all(engine)