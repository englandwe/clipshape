#!/usr/bin/env python

#Use this when you already have a mapped data file in clipmap format, but no clip data
#presumably you also have a txranges file, too.

#clip
#chr14_GL000194v1_random	90976	91001	RBM22_K562_rep01	200	-	
#chr start stop blah blah str

#map
#tx chr str pos value value

#USAGE: map clip

#TODO: stop hardcoding the chromosome dict
#General cleaning

import sys
import clipshape_core as csc
import copy

#functions

#translate clip chr names to shape/gft; catch mappings to chrs with no clip data
def fromClip(chr_name,chr_dict):
    if chr_name in chr_dict.keys():
        return chr_dict[chr_name]
    else:
        return 'noclip'

def buildMapDict(mapfile):
    mapdict={}
    with open(mapfile) as infile:
        for line in infile:
            tmpline=line.strip().split('\t')
            #chr,str,pos
            key=(tmpline[1],tmpline[2],tmpline[3])
            #assume nothing is a peak until proven otherwise
            try:
                mapdict[key].append(tmpline+[0])
            except KeyError:
                mapdict[key]=[tmpline+[0]]                
    return mapdict


##clip
#chr14_GL000194v1_random        90976   91001   RBM22_K562_rep01        200     -       
#chr start stop blah blah str
def addClip(clipfile,mapdict,chrdict):
    new_mapdict=copy.deepcopy(mapdict)
    with open(clipfile) as infile:
        for line in infile:
            tmpline=line.strip().split('\t')
            shapechr=fromClip(tmpline[0],chrdict)
            if shapechr != 'noclip':
                for pos in range(int(tmpline[1]),int(tmpline[2])):
                    key=(shapechr,tmpline[5],str(pos))
                    try:
                        for tx in new_mapdict[key]:
                            tx[-1] = 1
                    except KeyError:
                        sys.stderr.write('WARNING: %s not in mapfile\n' % str(key))
            else:
                sys.stderr.write('WARNING: %s not in chromosome dictionary\n' % tmpline[0])       
    return new_mapdict

#Chromosome tln dict 
#backwards in this case, so let's flip it       
chrnames={
    '1' : 'chr1',
    '2' : 'chr2',
    '3' : 'chr3',
    '4' : 'chr4',
    '5' : 'chr5',
    '6' : 'chr6',
    '7' : 'chr7',
    '8' : 'chr8',
    '9' : 'chr9',
    '10' : 'chr10',
    '11' : 'chr11',
    '12' : 'chr12',
    '13' : 'chr13',
    '14' : 'chr14',
    '15' : 'chr15',
    '16' : 'chr16',
    '17' : 'chr17',
    '18' : 'chr18',
    '19' : 'chr19',
    '20' : 'chr20',
    '21' : 'chr21',
    '22' : 'chr22',
    'MT' : 'chrM',
    'X' : 'chrX',
    'Y' : 'chrY'
}

chrnames_flipped={}
for key,val in chrnames.iteritems():
    chrnames_flipped[val] = key

#BEGIN

mapfile=sys.argv[1]
clipfile=sys.argv[2]

sys.stderr.write('building mapdict\n')
mapdict=buildMapDict(mapfile)

sys.stderr.write('commencing clipping\n')
new_mapdict=addClip(clipfile,mapdict,chrnames_flipped)

sys.stderr.write('writing output\n')
for vals in new_mapdict.itervalues():
    for val in vals:
        sys.stdout.write('\t'.join([str(x) for x in val])+'\n')
