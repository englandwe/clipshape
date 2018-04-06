#!/usr/bin/env python

#ENST00000425361        2       +       130190577       0.704   1       0       T       1       -18

import sys

motif=sys.argv[2]

subline=[]
with open(sys.argv[1]) as infile:
    prevline=next(infile).strip().split('\t')
    subline.append(prevline)
    for line in infile:
        tmpline=line.strip().split('\t')
        if tmpline[8] == prevline[8]:
            #consecutive
            subline.append(tmpline)
        else:
            #wrap up old range
            if subline[0][2] == '+':
                hasmotif=''.join([x[7].upper() for x in subline]).find(motif)
            elif subline[0][2] == '-':
                #subline=subline[::-1]
                #subline=[ x[0:9] + [-(int(x[9]))] for x in subline]
                hasmotif=''.join([x[7].upper() for x in subline]).find(motif)
            if hasmotif > -1:
                #set zero point
                offset=int(subline[hasmotif][9])
                for entry in subline:
                    entry.append(int(entry[9])-offset-2)
                    entry.append(motif)
                    print('\t'.join([str(x) for x in entry]))
            #new range
            subline=[]
            prevline=list(tmpline)
            subline.append(prevline)
#finish last range
subline.append(tmpline)
if subline[0][2] == '+':
    hasmotif=''.join([x[7].upper() for x in subline]).find(motif)
elif subline[0][2] == '-':
    #subline=subline[::-1]
    #subline=[x[0:9]+ [-(int(x[9]))] for x in subline]
    hasmotif=''.join([x[7].upper() for x in subline]).find(motif)
if hasmotif > -1:
    #set zero point
    offset=int(subline[hasmotif][9])
    for entry in subline:
        entry.append(int(entry[9])-offset-2)
        entry.append(motif)
        print('\t'.join([str(x) for x in entry]))
