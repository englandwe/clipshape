#!/usr/bin/env python

#USAGE: shapedata clipdata gtf output_prefix

#returns:  a file of transcript ranges (<output_prefix>.txranges)
#returns:  a file of mapped clip & shape data (<output_prefix>.clipmap)

#TODO: stop hardcoding the chromosome dict
#General cleaning


import sys
import clipshape_core as csc

#functions

#maps clip peak ranges to transcripts
def clipMapRanges(stitched_trans,clipdict):
    clipped_dict={}
    #only needed if clip and shape chrnames don't match
    clipchr=toClip(stitched_trans[1],chrnames)
    if clipchr == 'noclip':
        sys.stderr.write('WARNING: %s is mapped to %s, which has no CLIP data available. Skipping...\n' % (stitched_trans[0],stitched_trans[1]))
        return 'FAIL'
    else:
        for pos in stitched_trans[2]:
            isclip=0
            try:
                for rng in clipdict[clipchr]:
                    if pos >= rng[0] and pos <= rng[1]:
                        isclip=1
                        break
                if isclip==1:
                    clipped_dict[pos]=1
                else:
                    clipped_dict[pos]=0
            except KeyError:
                #no clip data for chr
                sys.stderr.write('WARNING: %s is mapped to %s, which has no CLIP data available. Skipping...\n' % (stitched_trans[0],stitched_trans[1]))
    return stitched_trans + [clipped_dict]

#grabs shape data for a stitched transcript.  takes output of clipMap
def clipMerge(clipped_trans,shape_list):
    for shape in shape_list:
        if shape[0] == clipped_trans[0]:
            shapevals=shape[3:]
            if len(shapevals) != len(clipped_trans[2]):
                sys.stderr.write('ERROR:  length mismatch in %s. Pos: %s Shape: %s.' % (shape[0],len(clipped_trans[2]),len(shapevals)) + os.linesep)
                return 'FAIL'
            else:
                return clipped_trans + [shapevals]
            break

#final merge
def mergeAll(shaped_trans):
    finallist=[]
    clipped_pos=shaped_trans[4].keys()
    for i in range(len(shaped_trans[2])):
        tmplist=[shaped_trans[0], shaped_trans[1], shaped_trans[3], shaped_trans[2][i], shaped_trans[5][i]]
        if tmplist[3] in clipped_pos:
            tmplist.append(shaped_trans[4][tmplist[3]])
        else:
            tmplist.append('0')
        finallist.append(tmplist)
    return finallist

#translate shape/gtf chr names to clip; catch mappings to chrs with no clip data
def toClip(chr_name,chr_dict):
    if chr_name in chr_dict.keys():
        return chr_dict[chr_name]
    else:
        return 'noclip'


###############################################
#Dictionary of shape/clip chr names
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

'''
chrnames={
    '1' : '1',
    '2' : '2',
    '3' : '3',
    '4' : '4',
    '5' : '5',
    '6' : '6',
    '7' : '7',
    '8' : '8',
    '9' : '9',
    '10' : '10',
    '11' : '11',
    '12' : '12',
    '13' : '13',
    '14' : '14',
    '15' : '15',
    '16' : '16',
    '17' : '17',
    '18' : '18',
    '19' : '19',
    '20' : '20',
    '21' : '21',
    '22' : '22',
    'MT' : 'MT',
    'X' : 'X',
    'Y' : 'Y'
}
'''
###############################################
#inputs

shapefile=sys.argv[1]
clipfile=sys.argv[2]
gtffile=sys.argv[3]
outprefix=sys.argv[4]

#import gtf
sys.stderr.write('INFO: importing gtf...\n')
gtfdict=csc.importGTF(gtffile)

#import clip data
sys.stderr.write('INFO: importing clip data...\n')
clip_data={}
with open(clipfile) as infile:
    for line in infile:
        cliptmp=line.strip().split('\t')
        try:
            clip_data[cliptmp[0]].append((int(cliptmp[1]),int(cliptmp[2])))
        except KeyError:
            clip_data[cliptmp[0]]=[(int(cliptmp[1]),int(cliptmp[2]))]

#import shape data
sys.stderr.write('INFO: importing SHAPE data...\n')
shape_data=[]
with open(shapefile) as infile:
    for line in infile:
        shape_data.append(line.strip().split('\t'))

############################################


ranges_outname = "%s.txranges" % outprefix
f1=open(ranges_outname,'w')

clip_outname = "%s.clipmap" % outprefix
f2=open(clip_outname,'w')

id_list=[x[0] for x in shape_data]
sys.stderr.write('INFO: Processing transcripts...\n')
txct=0
for id in id_list:
    trans_stitched=csc.stitchTrans(id,gtfdict)
    trans_ranges=csc.outputRanges(trans_stitched)
    trans_clipped=clipMapRanges(trans_stitched,clip_data)
    if trans_clipped != 'FAIL':
        trans_shaped=clipMerge(trans_clipped,shape_data)
        if trans_shaped != 'FAIL':
            trans_merged=mergeAll(trans_shaped)

            trans_out=csc.flattenList(trans_ranges)
            f1.write(trans_out+'\n')

            clip_out=csc.flattenList(trans_merged)
            f2.write(clip_out+'\n')
            txct+=1
            if txct % 100 == 0:
                sys.stderr.write('INFO: %s transcripts processed successfully...\n' % str(txct))

f1.close()
f2.close()
