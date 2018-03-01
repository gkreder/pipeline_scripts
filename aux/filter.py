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

def optimize_mz_tolerance(db):
	# Find max error among all transformations
	# look at all transformations between known observations
	# Starting with max error, incrementally decrease the mz tolerance until all 
	# the known-known transformations have been whitled out. 
	all_transformations = []
	known_transformations = []

	max_tol = float(db['transformations'].find_one()['dmz_err'])

	print('Optimizing MZ Tolerance...')

	for x in db['transformations'].all():
		all_transformations.append(x)
		
		if float(x['dmz_err']) > max_tol:
			max_tol = float(x['dmz_err'])

	best_tol = max_tol

	while 








def filter_transformations(args):
	db = dataset.connect('sqlite:///%s' % args.out_db)

	if os.path.exists(os.path.abspath(args.out_db)):
		print('Found output DB, overwriting...')
		os.system('rm \"%s\"' % os.path.abspath(args.out_db))

	os.system(('cp \"%s\" \"%s\"' % (os.path.abspath(args.in_db), (os.path.abspath(args.out_db)))))

	# edges = [x for x in db['edges'].all()]
	# transformations = {int(x['refNum']):x for x in db['transformations']}
	# refTransformations = {x['name']:x for x in db['refTransformations']}
	# knowns = {int(x['refNum']):x for x in db['knowns']}
	# observations = {int(x['refNum']):x for x in db['observations']}
	# nodes = {int(x['refNum']):x for x in db['nodes']}

	# print('Pruning impossible transformations')
	# impossible_count = 0
	# for trans in tqdm(transformations.values()):
	# 	obs_from = observations[trans['obs_from']]
	# 	obs_to = observations[trans['obs_to']]
	# 	refTrans = refTransformations[trans['trans']]
	# 	if ast.literal_eval(str(obs_from['known'])) and ast.literal_eval(str(obs_to['known'])):
	# 		if not isPossible(trans, obs_from, obs_to, refTrans):
	# 			# Remove Transformation
	# 			# Remove associated edges
	# 			db['transformations'].delete(refNum=trans['refNum'])
	# 			impossible_count += 1
	# 			for e in edges:
	# 				if int(e['trans']) == int(trans['refNum']):
	# 					db['edges'].delete(refNum=e['refNum'])
	# print('Removed %i of %i total transformations' % (impossible_count, len(transformations)))

	optimize_mz_tolerance(db)


