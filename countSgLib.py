#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gzip
import os
import sys
import datetime

def fileHandler(fileIn, operations = 'r', overwrite = False):
    if (os.path.exists(fileIn) and 'r'.upper() in operations.upper()) or ((not os.path.exists(fileIn) or (os.path.exists(fileIn) and overwrite)) and 'w'.upper() in operations.upper()):
        if fileIn.endswith('.gz'):
            try:
                fh = gzip.open(fileIn, operations.lower())
                return fh
            except Exception as e:
                print '\n[fileHandler]: error since [ %s ]! Now exit ...\n' % e
                sys.exit()
        else:
            try:
                fh = open(fileIn, operations.lower())
                return fh
            except Exception as e:
                print '\n[fileHandler]: error since [ %s ]! Now exit ...\n' % e
                sys.exit()
    else:
        print '\n[fileHandler]: could not proceed %s operation on %s since file not exists or not permitted to overwrite! Now exit ...\n' % (operations,fileIn)
        sys.exit()

def fastqReader(fh, numlines = 1):
    if fh.readable() or fh:
        return [ [ fh.readline().strip('\n') for x in xrange(4)] * numlines ]
    else:
        print '\n[fastqReader]: error, could not read file!\n'
        return [0]

def calcMatchAndMismatch(seq, primer, maxMis = 1):
    lenSeq = len(seq)
    lenPrimer = len(primer)
    pos = match = mis = 0
    for i in xrange(maxMis+1):
            pos = 0
            for j in xrange(lenSeq):
                    if lenSeq - lenPrimer - j + i < 0:
                            break
                    match = 0
                    mis = 0
                    for k in xrange(lenPrimer):
                            if i + k + 1 > lenPrimer:
                                    break
                            '''
                            Only allow primer overhang:
                              |--------------------------------    seq
                            -------------------->                  Primer
                            But not allow primer right overhang:
                            |---------------------------------     seq
                                              -------------------> Primer
                            The reson:
                            1. the in search-range sequence is long enough for primer identification;
                            2. the truncation in 5' of sequencing read is possible;
                            '''
                            if seq[j+k].upper() == primer[i+k].upper():
                                    match += 1
                            else:
                                    mis += 1
                                    if mis > maxMis:
                                            break
                    if mis <= maxMis and lenPrimer - mis == match:
                            pos = j - i
                            break
                    else:
                            continue
            if mis <= maxMis and lenPrimer - mis == match:
                    break
    if mis <= maxMis:
            return 1,pos,mis,match,pos+lenPrimer
    else:
            return 0

def matchWithMismatchPE(seq,primer,maxMis=1):
    pos1 = match1 = 0
    pos2 = match2 = 0
    seqArr = seq.split('-')
    seq1 = seqArr[0]
    seq2 = seqArr[1]
    primerArr = primer.split('-')
    primer1 = primerArr[0]
    primer2 = primerArr[1]
    lenPrimer1 = len(primer1)
    lenPrimer2 = len(primer2)
    if primer1 in seq1 and primer2 in seq2: # perfect match
            pos1 = seq1.index(primer1)
            match1 = lenPrimer1
            pos2 = seq2.index(primer2)
            match2 = lenPrimer2
            return 1,[match1,match2],[0,0],[pos1,pos2]
    match1 = calcMatchAndMismatch(seq1,primer1,maxMis)
    match2 = calcMatchAndMismatch(seq2,primer2,maxMis)
    if match1 and match2:
            if (match1[2] + match2[2]) <= maxMis*2 and match1[2] <= maxMis and match2[2] <= maxMis:
                    return 1,[match1[1],match2[1]],[match1[2],match2[2]],[match1[3],match2[3]],[match1[4],match2[4]]
    return 0

def matchWithMismatchSE(seq,primer,maxMis=1):
    match1 = 0
    seqArr = seq.split('-')
    seq1 = seqArr[0]
    primer1 = primer
    lenPrimer1 = len(primer)
    pos = ''
    if primer1 in seq1:
            pos = seq1.index(primer1)
            match = lenPrimer1
    if pos and pos != -1:
            return 1,[match],[0],[pos]
    match1 = calcMatchAndMismatch(seq1,primer1,maxMis)
    if match1:
            if match1[2] < maxMis:
                    return 1,[match1[1]],[match1[2]],[match1[3]],[match1[4]]
    return 0

def matchWithMismatch(seq,primer,maxMis=1):
    if "-" in primer:
            return matchWithMismatchPE(seq,primer,maxMis)
    else:
            return matchWithMismatchSE(seq,primer,maxMis)

def cacheKmer(inputType, klen = 7, byKey = False, removeFrequent = 10):
    uniqKmerDict = {}
    numUniq = 0
    if isinstance(inputType, dict):
        if byKey:
            for key in inputType.keys():
                for idx in range(len(key) - klen):
                    kmer = key[idx:idx+klen]
                    if not kmer in uniqKmerDict:
                        uniqKmerDict[kmer] = {'sample':[],'count':0}
                    uniqKmerDict[kmer]['sample'].append(inputType[key])
                    uniqKmerDict[kmer]['count'] += 1
        else:
            for key in inputType.keys():
                val = inputType[key]
                for idx in range(len(val) - klen):
                    kmer = val[idx:idx+klen]
                    if not kmer in uniqKmerDict:
                        uniqKmerDict[kmer] = {'sample':[],'count':0}
                    uniqKmerDict[kmer]['sample'].append(key)
                    uniqKmerDict[kmer]['count'] += 1
    elif isinstance(inputType, list):
        for seq in inputType:
            for idx in range(len(seq) - klen):
                    kmer = seq[idx:idx+klen]
                    if not kmer in uniqKmerDict:
                        uniqKmerDict[kmer] = {'sample':[],'count':0}
                    uniqKmerDict[kmer]['sample'].append(seq)
                    uniqKmerDict[kmer]['count'] += 1
    elif isinstance(inputType,str):
        for idx in range(len(inputType) - klen):
            kmer = inputType[idx:idx+klen]
            if not kmer in uniqKmerDict:
                uniqKmerDict[kmer] = {'sample':[],'count':0}
            uniqKmerDict[kmer]['sample'].append(inputType)
            uniqKmerDict[kmer]['count'] += 1
    else:
        print '\n[cacheKmer]: %s is not accepted for kmer cache!' % str(inputType)
        return 0
    for k,v in uniqKmerDict.items():
        if v['count'] == 1 or len(set(v['sample'])) == 1:
            numUniq += 1
        if len(set(v['sample'])) >= removeFrequent:
            uniqKmerDict.popitem(k)
        else:
            uniqKmerDict[k]['sample'] = list(set(v['sample']))
            uniqKmerDict[k]['count'] = len(set(v['sample']))
    return numUniq,uniqKmerDict

def fastDeter(kmerDictSeq,kmerDictLib):
    tmpSampleList = []
    for key in kmerDictSeq.keys():
        if key in kmerDictLib:
            tmpSampleList += kmerDictLib[key]['sample']
    return set(tmpSampleList)

def readSgFile(fh):
    sgDict = {}
    if fh:
        for line in fh:
            list1 = line.strip().split('\t')
            if len(list1) >= 2:
                if not list1[0] in sgDict:
                    sgDict[list1[0]] = list1[1]
                else:
                    print '\n[readSgFile]: %s - %s prior encountered for %s - %s!\n' % (list1[0],list1[1],list1[0],sgDict[list1[0]])
    return sgDict

def main():
    if len(sys.argv) != 4:
        print '\n[main]: need more parameters!\n\nTo run: python  %s  fastq1  fastq2  sgRNA_list\n' % sys.argv[0]
        sys.exit()
    fq1 = sys.argv[1]
    fq2 = sys.argv[2]
    sgFile = sys.argv[3]
    fh1 = fileHandler(fq1,'r')
    fh2 = fileHandler(fq2,'r')
    fh3 = fileHandler(sgFile,'r')
    sgDict1 = readSgFile(fh3)
    fh3.close()
    #sgList = sorted(sgDict1.keys(), key = lambda x:len(sgDict1[x]), reverse = True)
    matchCount = {}
    total_read = mapped = unmapped = 0
    memoryDict = {}
    numMemDict = 0
    sgRNAkmerDict = cacheKmer(sgDict1)
    while True:
        read1 = fastqReader(fh1,1)
        read2 = fastqReader(fh2,1)
        if not read1[0] or not read1[0][0]:
            print '\n[main]: read file finished!'
            break
        if total_read % 1000000 == 0:
            print '\n[%s] [main] read %.2f M reads ...' % (str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),1.0 * total_read / 1e6)
        total_read += 1
        matchList1 = []
        matchList2 = []
        status = 0
        for i in range(len(read1)):
            if '\t'.join((read1[i][1],read2[i][1])) in memoryDict:
                memoryDict['\t'.join((read1[i][1],read2[i][1]))]['count'] += 1
                if memoryDict['\t'.join((read1[i][1],read2[i][1]))]['mapped']:
                    mapped += 1
                else:
                    unmapped += 1
                matchCount[memoryDict['\t'.join((read1[i][1],read2[i][1]))]['key']] += 1
                read1.pop(i)
                read2.pop(i)
        for fqList1 in read1:
            read1seqkmerDict = cacheKmer(fqList1[1])
            sampleNameSet1 = fastDeter(read1seqkmerDict[1],sgRNAkmerDict[1])
            sampleNameList1 = sorted(sampleNameSet1, key = lambda x:len(sgDict1[x]), reverse = True)
            for sgName in sampleNameList1:
                match1 = matchWithMismatch(fqList1[1],sgDict1[sgName],1)
                if match1:
                    matchList1.append(sgName)
                    status += 1
                    break
            if status == 0:
                matchList1.append('')
        status = 0
        for fqList2 in read2:
            read2seqkmerDict = cacheKmer(fqList2[1])
            sampleNameSet2 = fastDeter(read2seqkmerDict[1],sgRNAkmerDict[1])
            sampleNameList2 = sorted(sampleNameSet2, key = lambda x:len(sgDict1[x]), reverse = True)
            for sgName in sampleNameList2:
                match1 = matchWithMismatch(fqList2[1],sgDict1[sgName],1)
                if match1:
                    matchList2.append(sgName)
                    status += 1
                    break
            if status == 0:
                matchList2.append('')
        for i in range(len(matchList1)):
            if not '\t'.join((read1[i][1],read2[i][1])) in memoryDict:
                memoryDict['\t'.join((read1[i][1],read2[i][1]))] = {'key':'','mapped':False,'count':0}
                numMemDict += 1
            if matchList1[i] and matchList2[i]: # both mapped
                mapped += 1
                memoryDict['\t'.join((read1[i][1],read2[i][1]))]['mapped'] = True
            else:
                unmapped += 1
            if '\t'.join((matchList1[i],matchList2[i])) in matchCount:
                matchCount['\t'.join((matchList1[i],matchList2[i]))] += 1
            else:
                matchCount['\t'.join((matchList1[i],matchList2[i]))] = 1
            memoryDict['\t'.join((read1[i][1],read2[i][1]))]['key'] = '\t'.join((matchList1[i],matchList2[i]))
            memoryDict['\t'.join((read1[i][1],read2[i][1]))]['count'] += 1
        if numMemDict >= 1000000:
            tmpKeyList = sorted(memoryDict.keys(), key = lambda x:memoryDict[x]['count'], reverse = True)
            medianCount = memoryDict[tmpKeyList[int(numMemDict*0.5)]]['count']
            afterDeleteKeys = { x:memoryDict[x] for x in tmpKeyList if memoryDict[x]['count'] > medianCount }
            memoryDict = afterDeleteKeys
            print '\n[main]: removed %d records with count less than %d from total %d stored keys ...' % (numMemDict - len(afterDeleteKeys.keys()),medianCount,numMemDict)
            numMemDict = len(afterDeleteKeys.keys())
    print '\n[%s] [main] read finished for %d reads ...' % (str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),total_read)
    fh1.close()
    fh2.close()
    print '\n#########################\n\nTotal reads:\t{}\nMapped:\t{}\t{}\nUnmapped:\t{}\t{}\n\n#########################\n'.format(total_read,mapped,('%.2f' % (100.0 * mapped / total_read)),unmapped,('%.2f' % (100.0 * unmapped / total_read)))
    for key in matchCount.keys():
        list2 = key.split('\t')
        tag = []
        for x in list2:
            if x:
                tag.append(x + '\t' + sgDict1[x])
            else:
                tag.append('-\t-')
        tag.append(str(matchCount[key]))
        print '\t'.join(tag)

if __name__ == '__main__':
    main()
