#!/usr/bin/env python

#ENST00000425361	2	+	130190577	0.704	1
#ENST00000425361	2	+	130190578	1.000	1

#ENST00000409276	2	-	38749704	0.000	1	53	G	11	-42	pos

import sys

pkfile=sys.argv[1]
localrange=range(-2,3)
#localrange=range(-7,8)

nazscores=[]
with open(pkfile) as infile:
    prevline=next(infile).strip().split('\t')
    startpos=prevline[3]
    #check if in desired local range
    if int(prevline[9]) in localrange:
        try:
            nazscores.append(float(prevline[4]))
        except ValueError:
            pass
    for line in infile:
        tmpline=line.strip().split('\t')
        if tmpline[0:3] == prevline[0:3] and int(tmpline[3]) in [int(prevline[3])+1,int(prevline[3])-1]:
            #consecutive
            if int(tmpline[9]) in localrange:
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
            if int(tmpline[9]) in localrange:
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
           


