#!/usr/bin/env python

#ENST00000425361        2       +       130190577       0.704   NEWVAL 1       0       T       1       -18
#motif should be a file with pwm in .pfm format (transposed from what comes out of homer)

import sys
import re
from Bio import motifs
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

def importMotif(motiffile):
    with open(motiffile) as infile:
        motif=motifs.read(infile,'pfm')
    return motif

motiffile=sys.argv[2]

motif=importMotif(motiffile)
motif_pssm=motif.pssm
#distribution = motif_pssm.distribution(precision=10**4)
#threshold = distribution.threshold_balanced(1000)
#sys.stderr.write("%s\n" % threshold)
threshold=7

subline=[]
with open(sys.argv[1]) as infile:
    prevline=next(infile).strip().split('\t')
    subline.append(prevline)
    for line in infile:
        tmpline=line.strip().split('\t')
        if tmpline[-2] == prevline[-2]:
            #consecutive
            subline.append(tmpline)
        else:
            #wrap up old range
            #if subline[0][2] == '+':
            #hasmotif=''.join([x[7].upper() for x in subline]).find(motif)
            #hasmotif=re.search(motif,''.join([x[7].upper() for x in subline]))
            hasmotif=[]
            #for position,score in motif_pssm.search(Seq(''.join([x[7].upper() for x in subline]),alphabet=IUPAC.unambiguous_dna),threshold=threshold):
            for position,score in motif_pssm.search(Seq(''.join([x[-3].upper() for x in subline]),alphabet=IUPAC.unambiguous_dna),threshold=threshold):
               #no reverse matches
               if position >= 0:
                   hasmotif.append([position,score])
            #take the best-scoring
            if len(hasmotif) > 0:
                best_hit=sorted(hasmotif, key=lambda x: x[1],reverse=True)[0]
                sys.stderr.write("%s,%s\n" % (best_hit[0],best_hit[1]))           
                #set zero point
                #offset=int(subline[best_hit[0]][9])
                offset=int(subline[best_hit[0]][-1])
                for entry in subline:
                    #entry.append(int(entry[9])-offset-2)
                    entry.append(int(entry[-1])-offset-2)
                    entry.append(motif.consensus)
                    print('\t'.join([str(x) for x in entry]))
            #new range
            subline=[]
            prevline=list(tmpline)
            subline.append(prevline)
#finish last range
subline.append(tmpline)
hasmotif=[]
#for position,score in motif_pssm.search(Seq(''.join([x[7].upper() for x in subline]),alphabet=IUPAC.unambiguous_dna),threshold=threshold):
for position,score in motif_pssm.search(Seq(''.join([x[-3].upper() for x in subline]),alphabet=IUPAC.unambiguous_dna),threshold=threshold):
    hasmotif.append([position,score])
    #take the best-scoring
if len(hasmotif) > 0:
    best_hit=sorted(hasmotif, key=lambda x: x[1],reverse=True)[0]
    #set zero point
    #offset=int(subline[best_hit[0]][9])
    offset=int(subline[best_hit[0]][-1])
    for entry in subline:
        #entry.append(int(entry[9])-offset-2)
        entry.append(int(entry[-1])-offset-2)
        entry.append(motif.consensus)
        print('\t'.join([str(x) for x in entry]))
