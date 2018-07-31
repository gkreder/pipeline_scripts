################################################################################
# Gabe Reder - gkreder@gmail.com
# 03/2018
# 
# Filtering a merged Sql database file
################################################################################
import sys
import os
import dataset
import aux.lib as lib
from tqdm import tqdm
################################################################################
# Returns a sorted list of observations, sorted by priority to look at
def rank(args):
	db = dataset.connect('sqlite:///%s' % args.in_db)
	db['observations'].create_index(['known'])
	db['transformations'].create_index(['obs_from', 'obs_to'])
	observations = [x for x in db['observations'].all()]
	knowns = [x for x in db.query('select * from observations where known')]
	unknowns = [x for x in db.query('select * from observations where not known')]
	transformations = [x for x in db['transformations'].all()]
	error_scores = {}
	guess_list = []
	for u_count, u in enumerate(tqdm(unknowns)):
		guesses = []
		guesses_transformations = {}
		u_num = u['refNum']
		name = u['name']
		query = "select * from transformations where obs_from=\'%s\' or obs_to=\'%s\'" % (u_num, u_num)
		# if name == 'u215_Taurine conjugation_75810':
			# print(query)
		K = [x for x in db.query(query)]
		S = 0
		if len(K) == 0:
			continue
		for t in K:
			S += t['dmz_err']
			if t['obs_from'] == u_num:
				k_num = t['obs_to']
			else:
				k_num = t['obs_from']
			k = db['observations'].find_one(refNum=k_num)
			refT = db['refTransformations'].find_one(name=t['trans'])
			if not k['known']:
				continue
			guess = lib.vecGuess(u, k, t, refT)
			guesses.append(guess)
			if lib.vecToString(guess, format = args.vec_format) in guesses_transformations:
				guesses_transformations[lib.vecToString(guess, format = args.vec_format)].append(t['refNum'])
			else:
				guesses_transformations[lib.vecToString(guess, format = args.vec_format)] = [t['refNum']]
			
		# if u_count > 2: 
		# 		break
		S /= len(K)
		S = 1 / S
		error_scores[name] = S
		guess_list.append((u, guesses, guesses_transformations))
	guess_list.sort(key = lambda tup : error_scores[tup[0]['name']], reverse = True)

	n_factor = 1.0 / sum(error_scores.values())
	for e in error_scores:
		error_scores[e] = (error_scores[e]*n_factor) + 1

	return guess_list, error_scores

def agreement_score(guesses):
	counts = []
	best_guess = guesses[0]
	best_count = 1
	for g in guesses:
		found = False
		for i, (d,count) in enumerate(counts):
			if g == d:
				counts[i] = (d, count + 1)
				found = True
				if count + 1 > best_count:
					best_guess = g
					best_count = count + 1
		if not found:
			counts.append((g, 1))

	# score = best_count - (len(counts) - 1) + 1
	score = best_count
	return score, best_guess


def clean(guess_list, error_scores, args):
	final_list = []
	agreement_scores = {}
	total_scores = {}

	for u, guesses, guesses_transformations in guess_list:
		if args.verbose:
			print('\n')
			print(u['name'])
			for g in guesses:
				print(g)
			# print(guesses)
			print('\n')

		if len(guesses) == 0:
			continue 

		a_score, best_guess = agreement_score(guesses)
		e_score = error_scores[u['name']]
		t_score = e_score**a_score
		# print(a_score, e_score, t_score)
		total_scores[u['name']] = t_score
		final_list.append((u, guesses, best_guess, guesses_transformations))

	final_list.sort(key = lambda tup : total_scores[tup[0]['name']], reverse = True)
	return final_list

def clean_rank(args):
	guess_list, error_scores = rank(args)
	if args.verbose:
		print('---------------------------------------')
		print('%i guess_list' % len(guess_list))
		print('---------------------------------------')
	final_list = clean(guess_list, error_scores, args)
	if args.verbose:
		print('---------------------------------------')
		print('%i final_list' % len(final_list))
		print('---------------------------------------')
		for r, (u, guesses, best_guess, guesses_transformations) in enumerate(final_list):
			print(u['name'], '(%i/%i)' % (r+1, len(final_list)) ,' --- ',lib.vecToString(best_guess))
	
	print('Writing results...')
	os.system('cp %s %s' % (args.in_db, args.out_db))
	out_db = dataset.connect('sqlite:///%s' % args.out_db)
	out_rows = []
	for r, (u, guesses, best_guess, guesses_transformations) in tqdm(enumerate(final_list)):
		out_row = {x: u[x] for x in u if x in ['name', 'mz', 'rt', 'refNum']}
		out_row['best_guess'] = lib.vecToString(best_guess, format = args.vec_format)
		out_row['guesses'] = str([lib.vecToString(g, format = 'pubchem') for g in guesses])
		out_row['rank'] = r + 1
		out_row['best_guess_transformations'] = ';'.join([x for x in guesses_transformations[out_row['best_guess']]])
		out_rows.append(out_row)

	out_db['unknowns'].insert_many(out_rows)

	print('Done')
	return final_list

# def write_cyto(args, final_list):
# 	print('Writing cytoscape file...')
# 	db = dataset.connect('sqlite:///%s' % args.in_db)
# 	transformations = [x for x in db.query('select * from transformations')]

# 	for r, (u, guesses, best_guess) in enumerate(tqdm(final_list)):
# 	# for r, (u, guesses, best_guess) in tqdm(enumerate(final_list)):
# 		# print(u, guesses, best_guess)
# 		print(u)
		
# 		if r > 2: 
# 			break
# 	# print(trans)
		







