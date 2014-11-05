#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import print_function
import getopt
import sys
import os
import stat
import re
import time
import platform
import subprocess

start_time = time.time()
exeDir = os.path.abspath(os.path.dirname(__file__))+os.sep

from itertools import izip, count
from Bio import SeqIO
from Bio import AlignIO

if os.path.exists(exeDir+"/tmp"):
	os.chmod(exeDir+"/tmp", stat.S_IRWXU+stat.S_IRWXG+stat.S_IRWXO)
	pass
else:
	os.mkdir(exeDir+"/tmp")
	os.chmod(exeDir+"/tmp", stat.S_IRWXU+stat.S_IRWXG+stat.S_IRWXO)

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
	elif species=="horse":
		for item in variantList:
			m=re.search(r'^(\d+)',item)
			if m:
				location=int(m.group(1))
				if location<16126 or location>16352:
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

def manualModifyVariant(species, variantList):
	variantNeedModify={"cattle":[["215+TCC"],["16200","16200+A"],["12171-12173d"]], "chicken":[["3941d"],["3940-3941d"]], "dog":[["15911d"],["796+T"],["16019+TGTAGCTGGAC"]], "goat":[["177G","179","180-181d"]], "horse":[["356d"],["5210d"],["16348-16355d"],["16403-16404d"],["355+C"],["16158-16357d"],["16334-16357d"]], "pig":[["136+C"],["136+CA"]], "sheep":[["16343","16343+C"],["567G","568C"],["15092G","15093C","15095A","15096C"]], "yak":[["880+GC"],["744+TC"],["893+GTGGGG"]]}
	variantAfterModify={"cattle":[["216","221+CCC"],["16199+A","16200"],["12173-12175d"]], "chicken":[["3941","3946d"],["3941","3945-3946d"]], "dog":[["15911","15918d"],["797","799+C"],["16019"]], "goat":[["177-178d","180T","181A"]], "horse":[["356","358d"],["5210","5217d"],["16127-16134d"],["16403d","16404d"],["356","356+T"],["16153-16352d"],["16329-16352d"]], "pig":[["137C","142+XA"],["137C","142+XA"]], "sheep":[["16342+C","16343"],["566+G","569d"],["15091+G","15097d"]], "yak":[["881G","892+CC"],["745","751+CC"],["892+GGTGGG"]]}
	beforeList=variantNeedModify[species]
	afterList=variantAfterModify[species]
	for index1, item1 in enumerate(beforeList):
		if len([item2 for item2 in item1 if item2 in variantList])==len(item1):
			for item in item1:
				variantList.remove(item)
			for item in afterList[index1]:
				variantList.append(item)
	return variantList

def runClustalwAndGetVariantSet(exeDir, species, clustalw):
	#begin to run clustalw
	p=subprocess.Popen('"'+exeDir+"/bin/"+clustalw+'"'+" -INFILE="+'"'+exeDir+"/tmp/in"+species+".fasta"+'"'+" -OUTFILE="+'"'+exeDir+"/tmp/out"+species+".fasta"+'"'+" -QUIET", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) 
	p.wait()
	#os.system(exeDir+"/bin/"+clustalw+" -INFILE="+exeDir+"/tmp/in"+species+".fasta"+" -OUTFILE="+exeDir+"/tmp/out"+species+".fasta"+" -QUIET")
	#analyse the output of clustalw
	tmpFile_handle=open(exeDir+"/tmp/out"+species+".fasta",'r')
	align = AlignIO.read(tmpFile_handle, "clustal")
	tmpFile_handle.close()
	refPrewChar=''
	queryPrewChar=''
	#parameter for insertation in query
	beginInsert = []
	endInsert = []
	insertString = []
	insertList = []
	totalInsertLong = 0
	#parameter for Ns in query
	beginNs = []
	endNs = []
	#parameter for deletion in query
	beginDel = []
	endDel = []
	delString = []
	#record variants
	variantList = []

	for num, refChar, queryChar in izip(count(),align[0].seq, align[1].seq):
		# query has inserts
		if refChar == '-':
			totalInsertLong = totalInsertLong + 1
			if refPrewChar != '-':
				beginInsert.append(num - totalInsertLong + 1)
				endInsert.append(num - totalInsertLong + 1)
				insertString.append(queryChar)
			else:
				endInsert[-1] = num - totalInsertLong + 1 
				insertString[-1] = insertString[-1] + queryChar
		if refChar != '-' and refChar != queryChar:
			# query has SNPs
			if queryChar != '-' and queryChar != 'N':
				if refChar == 'A':
					if queryChar == 'G':
						variantList.append(str(num - totalInsertLong+1))
					else:
						variantList.append(str(num - totalInsertLong+1)+queryChar)
				elif refChar == 'G':
					if queryChar == 'A':
						variantList.append(str(num - totalInsertLong+1))
					else:
						variantList.append(str(num - totalInsertLong+1)+queryChar)
				elif refChar == 'C':
					if queryChar == 'T':
						variantList.append(str(num - totalInsertLong+1))
					else:
						variantList.append(str(num - totalInsertLong+1)+queryChar)
				elif refChar == 'T':
					if queryChar == 'C':
						variantList.append(str(num - totalInsertLong+1))
					else:
						variantList.append(str(num - totalInsertLong+1)+queryChar)
				else:
					variantList.append(str(num - totalInsertLong+1)+queryChar)
			# query has Ns
			elif queryChar == 'N':
				if queryPrewChar != 'N':
					beginNs.append(num - totalInsertLong + 1)
					endNs.append(num - totalInsertLong + 1)
				else:
					endNs[-1] = num - totalInsertLong + 1
			# query has deletions
			else:
				if queryPrewChar != '-':
					beginDel.append(num - totalInsertLong + 1)
					endDel.append(num - totalInsertLong + 1)
					delString.append(refChar)
				else:
					endDel[-1] = num - totalInsertLong + 1
					delString[-1] = delString[-1] + refChar
		refPrewChar=refChar
		queryPrewChar=queryChar
	
	# the variants including SNPs, Ns and indels would be stored in querySet 
	querySet = set()
	# for insertion
	for beginPosition, string, endPosition in izip(beginInsert,insertString,endInsert):       
		querySet.add(str(beginPosition)+'+'+string)
	# for Ns
	for beginPosition, endPosition in izip(beginNs,endNs):
		if beginPosition == endPosition:
			querySet.add(str(beginPosition)+'N')
		else:
			querySet.add(str(beginPosition)+'-'+str(endPosition)+'N')
	# for deletion
	for beginPosition, string, endPosition in izip(beginDel,delString,endDel):
		if beginPosition == endPosition:
			querySet.add(str(beginPosition)+'d')
		else:
			querySet.add(str(beginPosition)+'-'+str(endPosition)+'d')
	querySet.update(variantList)
	return querySet

def analyzeSeqInFile(species, inputFile, output, exeDir, region):
	# read reference sequence and reference tree
	refFile = open(exeDir+"/ref/refSeq/"+species+"Ref.fasta", "r")
	refRecord = list(SeqIO.parse(refFile, "fasta"))[0]
	refFile.close()
	
	clustalw="clustalw2.exe"
	if platform.system() == "Linux":
		clustalw="clustalw2_linux"
	elif platform.system() == "Darwin":
		clustalw="clustalw2_mac"
	else:
		pass
	os.chmod(exeDir+"/bin/"+clustalw, stat.S_IRWXU+stat.S_IRWXG+stat.S_IRWXO)
	
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
	queryRecord = SeqIO.parse(queryFile, "fasta")
	queryNumber=0
	
	if output != "stdout":
		outFile=open(output, "w")
		print("Sample_Name\tHaplogroup\tMissing_Variants\tPrivate_Variants\tVariants", end="\n", file=outFile)
	else:
		print("Sample_Name\tHaplogroup\tMissing_Variants\tPrivate_Variants\tVariants", end="\n")

	# set parameter for modifying query sequence so as to give it correct sequence order consistent with refenece sequence
	cutRequireOne=False
	cutRequireTwo=False
	cutSeq=""
	#pig goat yak with wrong begin region (cut type one)
	if species=="pig":
		cutSeq="CAACCAAAACAAGCA"
		cutRequireOne=True
	elif species=="goat":
		cutSeq="GTTGATGTAGCTTAAA"
		cutRequireOne=True
	elif species=="yak":
		cutSeq="AACGCTATTAATATA"
		cutRequireOne=True
	#cattle with wrong end region (cut type two)
	elif species=="cattle":
		cutSeq="AAGACATCTCGATGG"
		cutRequireTwo=True
	else:
		pass
	
	for query in queryRecord:
		queryNumber=queryNumber+1
		output_handle = open(exeDir+"/tmp/in"+species+".fasta", "w")
		print(">refSeq_"+species+"\n"+refRecord.seq,end="\n", file=output_handle)
		
	#begin the modification of the order of query sequence
		#pig goat yak only
		R1=""
		R2=""
		if cutRequireOne and not query.seq.startswith(cutSeq):
			tmpos=query.seq.find(cutSeq)
			if tmpos>0:
				p1=query.seq[tmpos:]
				p2=query.seq[:tmpos]
				query.seq=p1+p2
		#cattle only
		elif cutRequireTwo and (not query.seq.endswith(cutSeq)) and region != "dloop":
			tmpos=query.seq.find(cutSeq)
			if tmpos>0:
				p1=query.seq[(tmpos+len(cutSeq)):]
				p2=query.seq[:(tmpos+len(cutSeq))]
				query.seq=p1+p2
		elif cutRequireTwo and region == "dloop":
			tmpos=query.seq.find(cutSeq)
			if tmpos>0:
				R1=query.seq[(tmpos+len(cutSeq)):]
				R2=query.seq[:(tmpos+len(cutSeq))]
		else:
			pass
		#end
		
		
		querySet=set()
		if cutRequireTwo and region == "dloop":
			#cattle dloop part1
			query.seq=R1
			print(">querySeq_"+species+"\n"+query.seq, end="\n", file=output_handle)
			output_handle.close()
			querySet1=runClustalwAndGetVariantSet(exeDir, species, clustalw)
			#cattle dloop part2
			output_handle = open(exeDir+"/tmp/in"+species+".fasta", "w")
			print(">refSeq_"+species+"\n"+refRecord.seq, end="\n", file=output_handle)
			query.seq=R2
			print(">querySeq_"+species+"\n"+query.seq, end="\n", file=output_handle)
			output_handle.close()
			querySet2=runClustalwAndGetVariantSet(exeDir, species, clustalw)
			querySet.update(querySet1)
			querySet.update(querySet2)
		else:
			print(">querySeq_"+species+"\n"+query.seq, end="\n", file=output_handle)
			output_handle.close()
			querySet=runClustalwAndGetVariantSet(exeDir, species, clustalw)

		# preprocessing querySet for pig and dog currently
		querySetTmp=set(removeVariantForNonSelectedRegion(removeVariantForSomeRegion(species,manualModifyVariant(species, list(querySet))), regionBegin, regionEnd, species, region))
		querySet={}
		querySet=querySetTmp.copy()
		querySetTmp={}
		# begin haplogroup assignment
		scoreMax = -100
		resultList = list() # store the assigned haplogroup
		missing = dict() # store the missing variants of each assigned haplogroup
		private = dict() # store the private variants of each assigned haplogroup
		for  index, haplogroupName in enumerate(treeHaplogroupNameList):
			refSet=treeHaplogroupVariantsSetList[index]
			score = 2*len(refSet & querySet) - len(refSet)
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
			print(query.id, end="\t", file=outFile)
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
			print(",".join(reorder(list(querySet))), end="\n", file=outFile)
		else:
			print(query.id, end="\t")
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
	-i or --input: Input file (a fasta file)
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
		region = re.sub(spaceReg, '', region) #/////
		#print(region)
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

