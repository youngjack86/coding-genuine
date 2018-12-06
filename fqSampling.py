#!/usr/bin/python
#coding: utf8

'''
Author: zhangyongjian@novogene.com
Date: Dec 6 2018
Restricted to personal use only

'''

import os
import sys
import gzip
import random

def fileHandler(file,operation='r'):
	if operation == 'r' and not os.path.isfile(file):
		print '[fileHandler]: the input file %s is not exist!' % file
		sys.exit()
	if operation == 'w' and os.path.isfile(file):
		print '[fileHandler]: the output file %s is exist!' % file
		sys.exit()
	if file.endswith('.gz'):
		fh = gzip.open(file,operation)
	else:
		fh = open(file,operation)
	return fh

def fqReader(fh,lines=1):
	fqList = []
	i = 0
	while i < 4 * lines:
		fqList.append(fh.readline().strip())
		i += 1
	return fqList

def getPrefix(file):
	basename = file.split('/')[-1]
	prefix = basename.split('.')[0]
	return prefix

def sampleStratgy(sampleRate=0.05, randSeed=1927):
	sampleBase = 10 ** len(str(sampleRate).split('.')[-1])
	sampleNumber = float(sampleRate) * sampleBase
	random.seed(randSeed)
	return([sampleBase,sampleNumber])

def randDetermine(sampleBase=100,sampleNumber=5):
	randNum = random.randint(1,sampleBase)
	if randNum <= sampleNumber:
		status = 1
	else:
		status = 0
	return status

def main():
	args = sys.argv
	if len(args) != 3:
		print '\n[main]: Need one parameter!\n\nTo RUN: python  %s  < fastq.gz >  < sampleRate >\n' % args[0]
		sys.exit()
	input = args[1]
	rate = args[2]
	output = 'sampled_' + getPrefix(input) + '.fq.gz'
	fhIN = fileHandler(input,'r')
	fhOUT = fileHandler(output,'w')
	randomSeed = 1927
	spstr = sampleStratgy(rate,randomSeed)
	spBase = spstr[0]
	spNum = spstr[1]
	totalSampledLength = 0
	lineNum = 0
	sampledFqLines = []
	maxCycle = 0
	numSampledLines = 0
	while True:
		fq = fqReader(fhIN)
		if fq[1] == '':
			if input.endswith('.gz'):
				fhIN.close()
				fhIN = fileHandler(input,'r')
			else:
				fhIN.seek(0,0)
			lineNum = 0
			maxCycle += 1
			if maxCycle > 10:
				break
			print '[main]: starts a new cycle of sampling [ cycle %d ], %d lines were sampled from last cycle!' % (maxCycle, numSampledLines)
			numSampledLines = 0
			continue
		lineNum += 1
		if not lineNum in sampledFqLines:
			lenRead = len(fq[1])
			stat = randDetermine(spBase,spNum)
			if stat > 0:
				sampledFqLines.append(lineNum)
				numSampledLines += 1
				totalSampledLength += lenRead
				fhOUT.write('\n'.join(fq) + '\n')
				if int(totalSampledLength / 1E6) != int((totalSampledLength - lenRead) / 1E6):
					print '\n\t%s bp sequences sampled!\n' % str(totalSampledLength)
		else:
			#print '[main]: %d exist in prior record. Skip it!' % lineNum
			pass
	fhIN.close()
	fhOUT.close()
	print '[main]: complished sampling, %d bp were sampled!\n' % totalSampledLength

if __name__ == '__main__':
	main()
