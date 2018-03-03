################################################################################
# Pipeline v2.0
# Gabe Reder - gkreder@gmail.com
# 01/2018
# 
# Filtering a merged Sql database file
################################################################################
import sys
import os
import dataset
from aux.lib import atomVector
import yaml
import ast
from tqdm import tqdm
import numpy as np
################################################################################
def isPossible(trans, obs_from, obs_to, refTrans):
	proposed_atom_change = yaml.load(refTrans['atoms'])
	atoms_from = atomVector(obs_from['SMILES'])
	atoms_to = atomVector(obs_to['SMILES'])
	possible = True
	for atom, atom_count in proposed_atom_change.items():
		if atom_count == 0:
			continue
		if atom_count < 0:
			if atom not in atoms_from.keys():
				possible = False
				continue
			if atoms_from[atom] < abs(atom_count):
				possible = False
				continue
		else:
			if atom not in atoms_to.keys() or atoms_to[atom] == 0:
				possible = False
				continue
			if atoms_to[atom] < atom_count:
				possible = False
				continue
	return possible

def no_id(x):
    del x['id']
    return x

def percentCorrect(known_transformations):
    correct = 0
    incorrect = 0
    for kt in known_transformations:
        if kt['known_correct']:
            correct += 1
        else:
            incorrect += 1
    return correct / (correct + incorrect)


def optimize_mz_tolerance(db):
	known_transformations = []
	edges = [x for x in db['edges'].all()]
	db['transformations'].create_index(['known_known', 'known_correct', 'dmz_err'])
	max_tol = float(db['transformations'].find_one()['dmz_err'])
	min_known_tol = float(db['transformations'].find_one()['dmz_err'])
	for x in db['transformations'].all():
	    if ast.literal_eval(str(x['known_known'])):
	        known_transformations.append(x)
	        if float(x['dmz_err']) < min_known_tol:
	            min_known_tol = float(x['dmz_err'])
	    if float(x['dmz_err']) > max_tol:
	        max_tol = float(x['dmz_err'])

	best_tol = max_tol
	best_percent = percentCorrect(known_transformations)
	t = max_tol
	increment = -1 * max_tol * 0.001
	temp_transformations = [x for x in known_transformations]

	print('Optimizing mz tolerance value (increment %f)...' % increment)
	# Find best tolerance value
	for t in tqdm(np.arange(max_tol, min_known_tol, increment)):
	    p = percentCorrect(temp_transformations)
	    if p > best_percent:
	        best_percent = p
	        best_tol = t
	    temp_transformations = [x for x in known_transformations if float(x['dmz_err']) < t]
	print('Done. Best tolerance = %f' % best_tol)

	delete_transformations = [no_id(x) for x in known_transformations if x['dmz_err'] > best_tol]
	keep_transformations = {x['refNum']:no_id(x) for x in known_transformations if x['dmz_err'] <= best_tol}
	keep_edges = [no_id(x) for x in edges if x['trans'] in keep_transformations.keys()]

	print('Deleting bad transformations...')
	db['transformations'].drop()
	db['edges'].drop()
	db.create_table('edges')
	db.create_table('transformations')
	db['transformations'].insert_many(keep_transformations.values())
	db['edges'].insert_many(keep_edges)




def filter_transformations(args):
	if os.path.exists(os.path.abspath(args.out_db)):
		print('Found output DB, overwriting...')
		os.system('rm \"%s\"' % os.path.abspath(args.out_db))

	os.system(('cp \"%s\" \"%s\"' % (os.path.abspath(args.in_db), (os.path.abspath(args.out_db)))))
	db = dataset.connect('sqlite:///%s' % args.out_db)

	edges = [x for x in db['edges'].all()]
	transformations = {x['refNum']:x for x in db['transformations']}
	refTransformations = {x['name']:x for x in db['refTransformations']}
	knowns = {x['refNum']:x for x in db['knowns']}
	observations = {x['refNum']:x for x in db['observations']}
	nodes = {x['refNum']:x for x in db['nodes']}

	print('Pruning impossible transformations')
	impossible_count = 0
	for trans in tqdm(transformations.values()):
		obs_from = observations[trans['obs_from']]
		obs_to = observations[trans['obs_to']]
		refTrans = refTransformations[trans['trans']]
		# if ast.literal_eval(str(obs_from['known'])) and ast.literal_eval(str(obs_to['known'])):
		if ast.literal_eval(str(trans['known_known'])):
			if not isPossible(trans, obs_from, obs_to, refTrans):
				# Remove Transformation
				# Remove associated edges
				db['transformations'].delete(refNum=trans['refNum'])
				impossible_count += 1
				for e in edges:
					if e['trans'] == trans['refNum']:
						db['edges'].delete(refNum=e['refNum'])
	print('Removed %i of %i total transformations' % (impossible_count, len(transformations)))

	optimize_mz_tolerance(db)

	# optimize_mz_tolerance(db)


