#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import print_function
import getopt
import sys
import os
import re
import time
import platform

start_time = time.time()
exeDir = os.path.abspath(os.path.dirname(__file__))+os.sep

from itertools import izip, count
from Bio import SeqIO
#from Bio import AlignIO
#from Bio.Align.Applications import ClustalwCommandline
#from Bio.Align.Applications import ClustalOmegaCommandline

def reorder(variantList):
	tmp=dict()
	for item in variantList:
		m=re.search(r'^(\d+)',item)
		if m:
			location=int(m.group(1))
			tmp.update({item:location})
	tmp2=sorted(tmp.items(),key = lambda tmp:tmp[1])
	variantListNew=[i[0] for i in tmp2]
	return variantListNew

def removeVariantForSomeRegion(species, variantList):
	variantListNew=list()
	if species=="dog":
		for item in variantList:
			m=re.search(r'^(\d+)',item)
			if m:
				location=int(m.group(1))
				if location<=16019:
					variantListNew.append(item)
		return variantListNew
	elif species=="pig":
		for item in variantList:
			m=re.search(r'^(\d+)',item)
			if m:
				location=int(m.group(1))
				if location<700 or location>921:
					variantListNew.append(item)
		return variantListNew
	else:
		return variantList

def removeVariantForNonSelectedRegion(variantList, regionBegin, regionEnd, species, region):
	variantListNew=list()
	if region =="dloop" and species=="cattle":
		for item in variantList:
			mD=re.search(r'^(\d+)-(\d+)d',item)
			m=re.search(r'^(\d+)',item)
			if mD:
				location1=int(mD.group(1))
				location2=int(mD.group(2))
				if location2 > regionBegin and location2 < regionEnd:
					variantListNew.append(str(location1)+"-"+str(regionBegin)+"d")
				elif location1 < regionEnd and location1 > regionBegin:
					variantListNew.append(str(regionEnd)+"-"+str(location2)+"d")
				elif location1 >= regionEnd:
					variantListNew.append(item)
				elif location2 <= regionBegin:
					variantListNew.append(item)
				else:
					pass
			elif m:
				location=int(m.group(1))
				if location<=regionBegin or location>=regionEnd:
					variantListNew.append(item)
			else:
				pass
		if "1-363d" in variantListNew:
			variantListNew.remove("1-363d")
		if "15792-16338d" in variantListNew:
			variantListNew.remove("15792-16338d")
		return variantListNew
	else:
		for item in variantList:
			m=re.search(r'^(\d+)',item)
			if m:
				location=int(m.group(1))
				if location>=regionBegin and location<=regionEnd:
					variantListNew.append(item)
		return variantListNew

def analyzeSeqInFile(species, inputFile, output, exeDir, region):
	# read reference sequence and reference tree
	refFile = open(exeDir+"/ref/refSeq/"+species+"Ref.fasta", "r")
	refRecord = list(SeqIO.parse(refFile, "fasta"))[0]
	refFile.close()
	
	regionSetDloop={"cattle":[363,15792],"chicken":[1,1232],"dog":[15461,16196],"horse":[15466,16657],"pig":[1,1254],"yak":[1,892],"sheep":[15437,16616],"goat":[15431,16642]} #the dloop of cattle is different from others
	regionSetWhole={"cattle":[1,16338],"chicken":[1,16785],"dog":[1,16196],"horse":[1,16657],"pig":[1,16690],"yak":[1,16322],"sheep":[1,16616],"goat":[1,16642]}
	regionSetNonDloop={"cattle":[364,15791],"chicken":[1233,16785],"dog":[1,15460],"horse":[1,15465],"pig":[1255,16690],"yak":[893,16322],"sheep":[1,15436],"goat":[1,15430]}
	regionBegin=0
	regionEnd=0
	treeFilePath = ""
	if region == "dloop":
		regionBegin=regionSetDloop[species][0]
		regionEnd=regionSetDloop[species][1]
		treeFilePath = exeDir+"/ref/treeFile/"+species+"Tree_dloop.txt"
	elif region == "nondloop":
		regionBegin=regionSetNonDloop[species][0]
		regionEnd=regionSetNonDloop[species][1]
		treeFilePath = exeDir+"/ref/treeFile/"+species+"Tree_withoutDloop.txt"
	elif region == "whole":
		regionBegin=regionSetWhole[species][0]
		regionEnd=regionSetWhole[species][1]
		treeFilePath = exeDir+"/ref/treeFile/"+species+"Tree_whole.txt"
	else:
		regionBeginTmp, regionEndTmp = region.split(":")
		if regionBeginTmp.isdigit() and regionEndTmp.isdigit():
			regionBegin=int(regionBeginTmp)
			regionEnd=int(regionEndTmp)
			if regionBegin < regionSetWhole[species][0] or regionEnd > regionSetWhole[species][1]:
				print("Region setting is out of the length of reference mtDNA, please check again.")
				sys.exit()
			elif regionBegin > regionEnd:
				print("Begin exceeds end for -r or --region parameter, please reset.")
				sys.exit()
			else:
				pass
		else:
			print("Region setting should use integers like 100:1000.")
			sys.exit()
		treeFilePath = exeDir+"/ref/treeFile/"+species+"Tree_whole.txt"

	treeFile = open(treeFilePath, "r")
	treeHaplogroupNameList=list()
	treeHaplogroupVariantsSetList=list(set())
	if region == "whole" or region =="dloop" or region =="nondloop":
		for line in treeFile:
			refHaplogrupName, refHaplogrupVariants = line.rstrip('\n').split(':')
			refVariantsSet=set(re.split(r',\s*',refHaplogrupVariants))
			treeHaplogroupNameList.append(refHaplogrupName)
			if len(refVariantsSet)!=0:
				treeHaplogroupVariantsSetList.append(refVariantsSet)
			else:
				treeHaplogroupVariantsSetList.append(set())
	else:
		for line in treeFile:
			refHaplogrupName, refHaplogrupVariants = line.rstrip('\n').split(':')
			refVariantsSet=set(removeVariantForNonSelectedRegion(re.split(r',\s*',refHaplogrupVariants), regionBegin, regionEnd, species, region))
			treeHaplogroupNameList.append(refHaplogrupName)
			if len(refVariantsSet)!=0:
				treeHaplogroupVariantsSetList.append(refVariantsSet)
			else:
				treeHaplogroupVariantsSetList.append(set())
	treeFile.close()

	queryFile = open(inputFile, "r")
	queryNumber=0
	
	if output != "stdout":
		outFile=open(output, "w")
		print("Sample_Name\tHaplogroup\tMissing_Variants\tPrivate_Variants\tMissing_of_your_judgement\tVariants", end="\n", file=outFile)
	else:
		print("Sample_Name\tHaplogroup\tMissing_Variants\tPrivate_Variants\tMissing_of_your_judgement\tVariants", end="\n")

	delFormat = re.compile(r'd\w*?$')
	spaceFormat = re.compile(r'\s+$')
	for query in queryFile:
		query=spaceFormat.sub('',query)
		query_id=""
		query_var=""
		query_userJudge=""
		inputList=list()
		inputList=re.split(r'\t+',query)
		if len(inputList)==2:
			query_id=inputList[0]
			query_var=inputList[1]
		elif len(inputList)==1:
			continue
		else:
			query_id=inputList[0]
			query_var=inputList[1]
			query_userJudge=inputList[2]
		# format input querySet
		queryList = list()
		querySet = set()
		for item in re.split(r',\s*',query_var):
			item = item.strip()
			item = item.replace('ins','+')
			item = delFormat.sub('d',item) # format deletion
			m = re.search(r'^(\D)(\d+)(\D)',item)
			if m:
				position = int(m.group(2)) - 1
				if position >= regionSetWhole[species][0] and position <= regionSetWhole[species][1]:
					refChar = refRecord.seq[position]
					queryChar = m.group(3)
					if refChar == 'A':
						if queryChar == 'G':
							item=m.group(2)
						else:
							item=m.group(2)+queryChar
					elif refChar == 'G':
						if queryChar == 'A':
							item=m.group(2)
						else:
							item=m.group(2)+queryChar
					elif refChar == 'C':
						if queryChar == 'T':
							item=m.group(2)
						else:
							item=m.group(2)+queryChar
					elif refChar == 'T':
						if queryChar == 'C':
							item=m.group(2)
						else:
							item=m.group(2)+queryChar
					else:
						item=m.group(2)+queryChar
			queryList.append(item)
		querySetTmp=set(removeVariantForNonSelectedRegion(removeVariantForSomeRegion(species,queryList), regionBegin, regionEnd, species, region))
		querySet=querySetTmp.copy()
		querySetTmp={}
		query_userJudgeList=re.split(r'\s*(,|:|;)\s*',query_userJudge)
		if len(querySet)==0:
			continue
		queryNumber=queryNumber+1
		# begin haplogroup assignment
		scoreMax = -100
		resultList = list() # store the assigned haplogroup
		missing = dict() # store the missing variants of each assigned haplogroup
		private = dict() # store the private variants of each assigned haplogroup
		missingUserJudgeOut = list() # store the missing variants of each haplogroup of user's assignment
		for  index, haplogroupName in enumerate(treeHaplogroupNameList):
			refSet=treeHaplogroupVariantsSetList[index]
			score = 2*len(refSet & querySet) - len(refSet)
			if haplogroupName in query_userJudgeList:
				tmpList= list(refSet - querySet)
				if len(tmpList)>0:
					missingUserJudgeOut.append(haplogroupName+": "+",".join(reorder(tmpList)))
				else:
					missingUserJudgeOut.append(haplogroupName+": "+"None")
			if score > scoreMax:
				scoreMax=score
				missing = {}
				private = {}
				resultList = []
				resultList.append(haplogroupName)
				privateTmp = list(querySet - refSet)
				missingTmp = list(refSet - querySet)
				if len(missingTmp) != 0:
					missing[haplogroupName] = missingTmp
				if len(privateTmp) != 0:
					private[haplogroupName] = privateTmp
			elif score == scoreMax:
				resultList.append(haplogroupName)
				privateTmp = list(querySet - refSet)
				missingTmp = list(refSet - querySet)
				if len(missingTmp) != 0:
					missing[haplogroupName] = missingTmp
				if len(privateTmp) != 0:
					private[haplogroupName] = privateTmp
			else:
				pass
		# output result
		if output != "stdout":
			print(query_id, end="\t", file=outFile)
			print(", ".join(resultList), end="\t", file=outFile)
			missingOut=list()
			privateOut=list()
			for haplogroup in resultList:
				if missing.get(haplogroup):
					missingOut.append(haplogroup+": "+",".join(reorder(missing[haplogroup])))
				if private.get(haplogroup):
					privateOut.append(haplogroup+": "+",".join(reorder(private[haplogroup])))
			print("; ".join(missingOut), end="\t", file=outFile)
			print("; ".join(privateOut), end="\t", file=outFile)
			print("; ".join(missingUserJudgeOut), end="\t", file=outFile)
			print(",".join(reorder(list(querySet))), end="\n", file=outFile)
		else:
			print(query_id, end="\t")
			print(", ".join(resultList), end="\t")
			missingOut=list()
			privateOut=list()
			for haplogroup in resultList:
				if missing.get(haplogroup):
					missingOut.append(haplogroup+": "+",".join(reorder(missing[haplogroup])))
				if private.get(haplogroup):
					privateOut.append(haplogroup+": "+",".join(reorder(private[haplogroup])))
			print("; ".join(missingOut), end="\t")
			print("; ".join(privateOut), end="\t")
			print("; ".join(missingUserJudgeOut), end="\t")
			print(",".join(reorder(list(querySet))), end="\n")
		#================
	print("--- Time cost: %d seconds ---" % (time.time() - start_time))
	print("--- Query number: %d ---" % queryNumber)
	queryFile.close()
	if output != "stdout":
		outFile.close()


try:
	options, args = getopt.getopt(sys.argv[1:],"vhs:r:i:o:",["version","help","species=","region=","input=","output="])
except getopt.GetoptError:
	sys.exit()

version="1.0" #current version
inputFile=""
inputDir=""
output="stdout" #default stdout
species=""
region="whole" #default whole


def printHelp():
	helpDoc='''\
MitoToolPy: a pipeline written in Python and designed for mitochondrial DNA analysis of domestic animal, for details, please see our website DomeTree (www.dometree.org)
Version: 1.0
License: GPLv3
Copyright: MitoTool & DomeTree Team
Contact: Long Fan (fan.long@mail.kiz.ac.cn)
Parameter:
	-v or --version: Software Version
	-h or --help: Help Document
	-s or --species: Setting a domestic animal from cattle, chicken, dog, goat, horse, pig, sheep, and yak.
	-r or --region: Setting region for analysis, e.g., whole, dloop, nondloop, 100:15000, (default: whole)
	-i or --input: Input file (a text file using tab characters as the delimiter of columns, the order of columns should be ID, Variants and Your_Haplogrouping(optional). Please see our demo for details.)
	-o or --output: Output file or stdout of console/terminal (default)'''
	print(helpDoc)
	
	
for name, value in options:
	if name in ("-h", "--help"):
		printHelp()
		sys.exit()
	elif name in ("-v", "--version"):
		print('Current version:',version)
		sys.exit()
	elif name in ("-s", "--species"):
		species=value.lower()
	elif name in ("-r", "--region"):
		region=value.lower()
		spaceReg = re.compile(r'\s+')
		region = re.sub(spaceReg, '', region)
	elif name in ("-i", "--input") and os.path.exists(value):
		if os.path.isfile(value):
			inputFile=value
		elif os.path.isdir(value):
			inputDir=value
		else:
			pass
	elif name in ("-o", "--output"):
		if not os.path.exists(value) and os.path.exists(os.path.split(value)[0]) and value!="stdout":
			output=value
		elif not os.path.exists(value) and not os.path.exists(os.path.split(value)[0]) and value!="stdout":
			print("The directory of output file does not exists, please check again.")
			sys.exit()
		elif os.path.isdir(value):
			print("Output file should be a file but not a directory, please check again.")
			sys.exit()
		elif os.path.exists(value) and os.path.isfile(value):
			checking=raw_input("Output file already exists, whether to overwrite it? (y/n)")
			checking=checking.lower()
			print()
			if checking in ("y","yes"):
				output=value
			else:
				sys.exit()
		else:
			pass
	else:
		pass

if inputFile =="" and inputDir=="":
	print("Please set an input file or an directory contains batch input files using -i or --input.")
	sys.exit()
		
if species not in ("cattle", "chicken", "dog", "goat", "horse", "pig", "sheep", "yak"):
	print("Please set an available domestic animal using -s or --species, the candidates include cattle, chicken, dog, goat, horse, pig, sheep, and yak.")
	sys.exit()
	
print("============ Parameter ============")
if inputFile!="":
	print('python "'+os.path.abspath(__file__)+'"', " -s "+species, " -r "+'"'+region+'"', " -i "+'"'+inputFile+'"', " -o "+'"'+output+'"\n')
else:
	print('python "'+os.path.abspath(__file__)+'"', " -s "+species, " -r "+'"'+region+'"', " -i "+'"'+inputDir+'"', " -o "+'"'+output+'"\n')
print("========= Program running =========")

# python C:\Users\user\Desktop\dometree-tool\mitotoolpy.py -s dog -i C:\Users\user\Desktop\dometree-tool\test\dog.fasta -o C:\Users\user\Desktop\dometree-tool\test\dog_test.txt -r "1:100"
analyzeSeqInFile(species, inputFile, output, exeDir, region)