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

	# Check correctness of edges
	edges = [x for x in db['edges'].all()]
	transformations = {int(x['refNum']):x for x in db['transformations']}
	refTransformations = {x['name']:x for x in db['refTransformations']}
	knowns = {int(x['refNum']):x for x in db['knowns']}
	observations = {int(x['refNum']):x for x in db['observations']}


	e_corrects = []
	e_incorrects = []
	skipped_edges = 0
	print('')
	print('Validating edges...')
	for e in tqdm(edges):
		trans = transformations[int(e['trans'])]
		refTrans = refTransformations[trans['trans']]
		obs_from = observations[int(trans['obs_from'])]
		obs_to = observations[int(trans['obs_to'])]
		if not (ast.literal_eval(str(trans['known_known']))):
			skipped_edges += 1
			continue

		form_from = obs_from['SMILES']
		form_to = obs_to['SMILES']
		v_from = atomVector(form_from)
		v_to = atomVector(form_to)
		v_diff = vecDifference(v_from, v_to)
		v_expected = yaml.load(refTrans['atoms'])
		atoms_observed = v_diff
		atoms_expected = v_expected
		if vecEquals(v_diff, v_expected):
			e_corrects.append({'edge' : e, 'transformation' : trans, 'observed_vector' : v_diff, 'expected_vector' : v_expected})
		else:
			e_incorrects.append({'edge' : e, 'transformation' : trans, 'observed_vector' : v_diff, 'expected_vector' : v_expected})

	# Check correctness of all transformations
	transformations = [x for x in db['transformations'].all()]

	t_corrects = []
	t_incorrects = []
	skipped_transformations = 0
	print('Validating transformations...')
	for trans in tqdm(transformations):
		refTrans = refTransformations[trans['trans']]
		obs_from = observations[int(trans['obs_from'])]
		obs_to = observations[int(trans['obs_to'])]
		if not (ast.literal_eval(str(trans['known_known']))):
			skipped_edges += 1
			continue
		form_from = obs_from['SMILES']
		form_to = obs_to['SMILES']
		v_from = atomVector(form_from)
		v_to = atomVector(form_to)
		v_diff = vecDifference(v_from, v_to)
		v_expected = yaml.load(refTrans['atoms'])
		atoms_observed = v_diff
		atoms_expected = v_expected
		if vecEquals(v_diff, v_expected):
			t_corrects.append({'transformation' : trans, 'observed_vector' : v_diff, 'expected_vector' : v_expected})
		else:
			t_incorrects.append({'transformation' : trans, 'observed_vector' : v_diff, 'expected_vector' : v_expected})

	e_correct_percent = len(e_corrects) / len(e_corrects + e_incorrects)
	t_correct_percent = len(t_corrects) / len(t_corrects + t_incorrects)

	print('')
	print(len(t_corrects), len(t_incorrects), skipped_transformations)
	print('%f correct transformations ' % t_correct_percent)
	print('')
	print(len(e_corrects), len(e_incorrects), skipped_edges)
	print('%f correct edges ' % e_correct_percent)

	results = {'edges_correct_percent' : e_correct_percent, 
	'transformations_correct_percent' : t_correct_percent,
	'correct_edges' : e_corrects,
	'incorrect_edges' : e_incorrects,
	'correct_transformations' : t_corrects,
	'incorrect_transformations' : t_incorrects}

	return results

