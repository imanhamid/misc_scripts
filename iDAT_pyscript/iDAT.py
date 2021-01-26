#! /usr/bin/env python

import msprime, pyslim
import json
import numpy as np
import pandas as pd
import math
import sys
import re
from sklearn.metrics import auc

def dumpTracts(samples, anc_start, anc_end):
	'''
	Returns EdgeTable mapping P1 or P2 ancestry across all chromosomes
	'''
	tablecoll = ts.tables
	ancestor_map = tablecoll.map_ancestors(samples, range(anc_start, anc_end))
	return ancestor_map

def getChildTracts(child, EdgeTable):
	'''
	Returns array of tract starts and ends for one chromosome
	'''
	indices = np.argwhere(EdgeTable.child==child)
	child_arr =[]
	for index in indices:
		child_arr.append([EdgeTable[int(index)].left, EdgeTable[int(index)].right])
		child_arr.sort(key = lambda x: x[0])
	return child_arr

def getMergedTracts(child_arr):
	'''
	Returns merged tracts for one chromosome
	'''
	tracts = []
	start = -1
	max = -1
	for i in range(len(child_arr)):
		t = child_arr[i]
		if t[0] > max:
			if i != 0:
				tracts.append([start, max])
			max = t[1]
			start = t[0]
		else:
			if t[1] >= max:
				max = t[1]
	if max != -1 and [start,max] not in tracts:
		tracts.append([start,max])
	return tracts

def beneTract(merged_tracts, locus):
	'''
	Returns tract length, start, and end for tract spanning a locus
	'''
	bene_tract_index = int(np.argwhere(np.logical_and(np.array(merged_tracts)[:,0] <= locus,
		locus <= np.array(merged_tracts)[:,1])))

	tract = merged_tracts[bene_tract_index]

	if tract[0] < chr_lengths[chr_index][0]:
		tract_start = chr_lengths[chr_index][0]
	else:
		tract_start = tract[0]
	if tract[1] > chr_lengths[chr_index][1]:
		tract_end = chr_lengths[chr_index][1]
	else:
		tract_end = tract[1]
	distance1 = int(tract_end - locus)
	distance2 = int(locus - tract_start)

	return (distance1, distance2)

def getiDAT(locus, sample_tracts):
	abs_dist = [abs(locus-i) for i in chr_coords]
	abs_dist.sort()
	abs_dist=np.array(abs_dist)
	DAT = np.array([(sum(i >= x for i in sample_tracts)/len(sample_tracts))**2 for x in abs_dist])
	try:
		iDAT = auc(abs_dist[DAT>=0.25], DAT[DAT>=0.25])
		return(iDAT)
	except ValueError:
		return("NA")

if __name__=="__main__":
	infile = sys.argv[1]
	locus = int(sys.argv[2])
	outfile = re.search("(.*).trees", infile).group(1)

	chr_lengths = np.array([[0,249904549],[249904550,493103922],
		[493103923,691126352],[691126353,882661886],[882661887,1063577146],
		[1063577147,1234692213],[1234692214,1394013772],[1394013773,1540453883],
		[1540453884,1682150456],[1682150457,1817685203],[1817685204,1952731822],
		[1952731823,2086583717],[2086583718,2201753595],[2201753596,2309103135],
		[2309103136,2411634527],[2411634528,2501989280],[2501989281,2583518887],
		[2583518888,2661600397],[2661600398,2720981238],[2720981239,2784006758],
		[2784006759,2832164335],[2832164336,2883468901]])

	chr_index = int(np.argwhere(np.logical_and(chr_lengths[:,0]<=locus, locus<=chr_lengths[:,1])))
	
	chrom = chr_index+1

	coords_file = f"/work/ih49/chr{chrom}_coords.json"

	with open(coords_file) as f:
		chr_coords = json.load(f)

	ts = pyslim.load(infile).simplify()

	sample_high = ts.num_samples

	rng = np.random.default_rng()
	sample1 = rng.choice(range(4, sample_high, 2), size=172, replace=False)
	sample2 = sample1 + 1
	roots_sample = np.array([0,1,2,3])
	sample3 = np.concatenate([roots_sample, sample1, sample2])

	ts = ts.simplify(sample3)

	bene_tree = ts.at(locus)

	#get samples with P1 ancestry at bene_locus
	bene_samples = []

	for root in bene_tree.roots:
		if bene_tree.population(root)==0:
			for sample in bene_tree.samples(root):
				if sample >= 4:
					bene_samples.append(sample)

	#get samples with P2 ancestry at bene_locus
	non_bene_samples = []

	for root in bene_tree.roots:
		if bene_tree.population(root)==1:
			for sample in bene_tree.samples(root):
				if sample >= 4:
					non_bene_samples.append(sample)

	if len(bene_samples) and len(non_bene_samples):
		iDAT = []
		for anc in ["P1", "P2"]:
			if anc == "P1":
				anc_start = 0
				anc_end = 2
				samples = bene_samples
				sample_tracts = []
			elif anc == "P2":
				anc_start = 2
				anc_end = 4
				samples = non_bene_samples
				sample_tracts = []
			anc_map = dumpTracts(samples, anc_start, anc_end)
			for child in samples:
				sample_tracts.extend(beneTract(getMergedTracts(getChildTracts(child, anc_map)), locus))
			iDAT.append(getiDAT(locus, sample_tracts))
		iDAT_score = np.log(iDAT[1]/iDAT[0])
		print(iDAT_score)
	else:
		print("NA")
