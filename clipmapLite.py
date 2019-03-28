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

#functions

#translate shape/gtf chr names to clip; catch mappings to chrs with no clip data
def toClip(chr_name,chr_dict):
    if chr_name in chr_dict.keys():
        return chr_dict[chr_name]
    else:
        return 'noclip'


##clip
#chr14_GL000194v1_random        90976   91001   RBM22_K562_rep01        200     -       
#chr start stop blah blah str
def buildClipDict(clipfile):
    clipdict={}
    with open(clipfile) as infile:
        for line in infile:
            tmpline=line.strip().split('\t')
            posrange=[(x,tmpline[5]) for x in range(int(tmpline[1]),int(tmpline[2]))]
            try:
                clipdict[tmpline[0]]+=posrange
            except KeyError:
                clipdict[tmpline[0]]=posrange
    return clipdict

def addClip(mapfile,clipdict,chrdict):
    linect=0
    with open(mapfile) as infile:
        for line in infile:
            linect+=1
            if linect % 1000 == 0:
                sys.stderr.write('Processing line %s...\n' % linect)
            tmpline=line.strip().split('\t')
            clipchr=toClip(tmpline[1],chrdict)
            if clipchr != 'noclip':
                if (tmpline[3],tmpline[2]) in clipdict[clipchr]:
                    isPeak=1
                else:
                    isPeak=0
                outline=tmpline+[isPeak]
                print('\t'.join([str(x) for x in outline]))
            else:
                sys.stderr.write('WARNING: %s not in clip dictionary\n' % tmpline[1])



#Chromosome tln dict        
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

#BEGIN

mapfile=sys.argv[1]
clipfile=sys.argv[2]

sys.stderr.write('building clipdict\n')
clipdict=buildClipDict(clipfile)

sys.stderr.write('commencing clipping\n')
addClip(mapfile,clipdict,chrnames)

