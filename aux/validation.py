################################################################################
# Gabe Reder - gkreder@gmail.com
# Validation aux script
################################################################################
import sys
import os
import dataset
import pickle as pkl
from aux.lib import Obs, Node, MergedSet, RefTrans, atomVector, vecDifference, vecEquals
import xlsxwriter
import yaml
from tqdm import tqdm
import ast
################################################################################
def validate(args):
	db = dataset.connect('sqlite:///%s' % args.db)

	results = {
	'correct_edges' : [],
	'incorrect_edges' : [],
	'skipped_edges' : [],
	'correct_transformations' : [],
	'incorrect_transformations' : [],
	'skipped_transformations' : []
	}


	# Check Edges
	edges = [x for x in db['edges'].all()]
	transformations = {x['refNum'] : x for x in db['transformations'].all()}
	# transformations = {}
	# for x in db['transformations'].all():
	# 	print(x['refNum'])
	# 	transformations[x['refNum']] = x


	num_correct_edges = 0
	num_incorrect_edges = 0
	skipped_edges = 0
	for e in edges:
		refNum = e['trans']
		trans = transformations[refNum]
		if not ast.literal_eval(str(trans['known_known'])):
			skipped_edges += 1
			results['skipped_edges'].append(e)
			continue
		if ast.literal_eval(str(trans['known_correct'])):
			num_correct_edges += 1
			results['correct_edges'].append(e)
		else:
			num_incorrect_edges += 1
			results['incorrect_edges'].append(e)

	edges_correct_percent = num_correct_edges / (num_correct_edges + num_incorrect_edges)
	results['edges_correct_percent'] = edges_correct_percent

	# Check Transformations
	transformations = [x for x in db['transformations'].all()]
	skipped_transformations = 0
	num_correct_transformations = 0
	num_incorrect_transformations = 0
	for trans in transformations:
		if not ast.literal_eval(str(trans['known_known'])):
			skipped_transformations += 1
			results['skipped_transformations'].append(trans)
			continue
		if ast.literal_eval(str(trans['known_correct'])):
			num_correct_transformations += 1
			results['correct_transformations'].append(trans)
		else:
			num_incorrect_transformations += 1
			results['incorrect_transformations'].append(trans)

	transformations_correct_percent = num_correct_transformations / (num_correct_transformations + num_incorrect_transformations)
	results['transformations_correct_percent'] = transformations_correct_percent

	print('Correct Transformations - %i' % num_correct_transformations)
	print('Incorrect Transformations - %i' % num_incorrect_transformations)
	print('Percent Correct Transformations %f' % transformations_correct_percent)
	print('')
	print('Correct Edges - %i' % num_correct_edges)
	print('Incorrect Edges - %i' % num_incorrect_edges)
	print('Percent Correct Edges %f' % edges_correct_percent)

	return(results)