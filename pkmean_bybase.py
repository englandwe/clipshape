#!/usr/bin/env python

#ENST00000318388	18	+	9122581	0.055	1	0	T	6	-47	-52	GCATG

import sys

desired_bases=sys.argv[2].strip().split(',')
#print(desired_bases)
#localrange=range(-2,3)
localrange=range(-7,8)

nazscores=[]
with open(sys.argv[1]) as infile:
    prevline=next(infile).strip().split('\t')
    startpos=prevline[3]
    #print(prevline[7])
    if int(prevline[9]) in localrange and prevline[7] in desired_bases:
        #print("in")
        try:
            nazscores.append(float(prevline[4]))
        except ValueError:
            pass
    for line in infile:
        tmpline=line.strip().split('\t')
        if tmpline[0:3] == prevline[0:3] and int(tmpline[3]) in [int(prevline[3])+1,int(prevline[3])-1]:
            #consecutive
            #print(tmpline[7])
            if int(tmpline[9]) in localrange and tmpline[7] in desired_bases:
                #print("in")
                try:
                    nazscores.append(float(tmpline[4]))
                except ValueError:
                    pass
            prevline=list(tmpline)
        else:
            #wrap up old range
            try:
                nazmean=sum(nazscores)/float(len(nazscores))
            except ZeroDivisionError:
                nazmean='NaN'
            stoppos=prevline[3]
            print('\t'.join(prevline[0:3]+[startpos,stoppos,str(nazmean),prevline[5],prevline[10],str(len(nazscores))]))
            #new range
            nazscores=[]
            startpos=tmpline[3]
            if int(tmpline[9]) in localrange and tmpline[7] in desired_bases:
                try:
                    nazscores.append(float(tmpline[4]))
                except ValueError:
                    pass
            prevline=list(tmpline)
#finish last range
    try:
        nazmean=sum(nazscores)/float(len(nazscores))
    except ZeroDivisionError:
        nazmean='NaN'
    stoppos=prevline[3]
    print('\t'.join(prevline[0:3]+[startpos,stoppos,str(nazmean),prevline[5],prevline[10],str(len(nazscores))]))
           



