#!/usr/bin/env python

#Core classes and functions for the clipshape pipeline

import sys
import os
import random
import re
from Bio import SeqIO,motifs
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from itertools import groupby
from operator import itemgetter

#classes

#handle gtfs with ease and grace
class GtfRec(object):
    def __init__(self,reclist):
         self.seqname=reclist[0]
         self.source=reclist[1]
         self.feature=reclist[2]
         self.start=int(reclist[3])
         self.end=int(reclist[4])
         self.score=reclist[5]
         self.strand=reclist[6]
         self.frame=reclist[7]
         self.attdict={}
         for self.item in reclist[8].strip(';').split('; '):
             self.splitline=self.item.replace('\"','').split(' ')
             self.attdict[self.splitline[0]]=self.splitline[1]

#functions

#understands gencode gtf files
#imports everything with transcript id (i.e. not a gene)
#puts it in a dict by transcript id
def importGTF(gtffile):
    gtfdict={}   
    skipcount=0    
    sys.stderr.write('INFO: Importing gtf: %s.\n' % gtffile)                
    with open(gtffile) as infile:
        for line in infile:
            if not line.startswith('#'):
                rec=GtfRec(line.strip().split('\t'))
                try:
                    gtfkey=rec.attdict['transcript_id']
                    try:
                        gtfdict[gtfkey].append(rec)
                    except KeyError:
                        gtfdict[gtfkey]=[rec]
                except KeyError:
                    #no transcript id == it was probably a gene
                    skipcount+=1
    sys.stderr.write('INFO: gtf imported; %s records had no transcript id and were skipped.\n' % str(skipcount))                
    return gtfdict

#stitch together exons to make a transcript
def stitchTrans(ensembl_id,gtfdict):
    pos_list=[]
    chrid=''
    tmplist=[]
    for gtf_rec in gtfdict[ensembl_id]:
        #grab each exon's start & stop, plus exon number so they can be put in order
        if gtf_rec.feature == 'exon':
            tmplist.append([gtf_rec.attdict['exon_number'],gtf_rec.start,gtf_rec.end+1])
        if len(chrid) == 0:
            chrid=gtf_rec.seqname
            strand=gtf_rec.strand
    #sort the list by exon number
    sorted_tmplist=sorted(tmplist, key=lambda x: int(x[0]))
    for item in sorted_tmplist:
        pos_list+=range(item[1],item[2])
    return [ensembl_id,chrid,pos_list,strand]

#takes a stitched transcript and spits out a set of ranges
def outputRanges(stitched_trans):
    subrange_list=[]
    outlist=[]
    for k, g in groupby(enumerate(stitched_trans[2]), lambda (i, x): i-x):
        subrange_list.append(map(itemgetter(1), g))
    for subrange in subrange_list:
        sub_start=min(subrange)
        sub_stop=max(subrange)
        outlist.append([stitched_trans[0],stitched_trans[1],stitched_trans[3],sub_start,sub_stop])
    return outlist

#exactly what it says
def flattenList(listin):
    list2=[]
    for item in listin:
        list2.append('\t'.join([str(x) for x in item]))
    final='\n'.join(list2)
    return final
