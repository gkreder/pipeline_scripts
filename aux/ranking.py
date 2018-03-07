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
	for u in tqdm(unknowns):
		guesses = []
		u_num = u['refNum']
		name = u['name']
		query = "select * from transformations where obs_from=\'%s\' or obs_to=\'%s\'" % (u_num, u_num)
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
		S /= len(K)
		S = 1 / S
		error_scores[name] = S
		guess_list.append((u, guesses))
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

	for u, guesses in guess_list:
		if len(guesses) == 0:
			continue 

		a_score, best_guess = agreement_score(guesses)
		e_score = error_scores[u['name']]
		t_score = e_score**a_score
		# print(a_score, e_score, t_score)
		total_scores[u['name']] = t_score
		final_list.append((u, guesses, best_guess))

	final_list.sort(key = lambda tup : total_scores[tup[0]['name']], reverse = True)
	return final_list

def clean_rank(args):
	guess_list, error_scores = rank(args)
	final_list = clean(guess_list, error_scores, args)
	for (u, guesses, best_guess) in final_list:
		print(u['name'], best_guess)
	return final_list







