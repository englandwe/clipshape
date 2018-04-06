#!/usr/bin/env python

#grab an eclip bedfile and extract the relevant bits

import sys
import clipshape_core as csc

chrnames={
 'chr1': '1',
 'chr10': '10',
 'chr11': '11',
 'chr12': '12',
 'chr13': '13',
 'chr14': '14',
 'chr15': '15',
 'chr16': '16',
 'chr17': '17',
 'chr18': '18',
 'chr19': '19',
 'chr2': '2',
 'chr20': '20',
 'chr21': '21',
 'chr22': '22',
 'chr3': '3',
 'chr4': '4',
 'chr5': '5',
 'chr6': '6',
 'chr7': '7',
 'chr8': '8',
 'chr9': '9',
 'chrM': 'MT',
 'chrX': 'X',
 'chrY': 'Y'
 }

def tlbed(bedfile,chrdict):
    outlist=[]
    with open(bedfile) as infile:
        for line in infile:
            tmpline=line.strip().split('\t')
            try:
                outline=[chrdict[tmpline[0]]]+tmpline[1:3]+[tmpline[5]]
                outlist.append(outline)
            except KeyError:
                sys.stderr.write('WARNING: %s not in chromosome dictionary.  Skipping...\n' % tmpline[0])
    return outlist

bedfile=sys.argv[1]
print(csc.flattenList(tlbed(bedfile,chrnames)))
