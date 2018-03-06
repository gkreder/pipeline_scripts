################################################################################
# Gabe Reder - gkreder@gmail.com
# 03/2018
# 
# Filtering a merged Sql database file
################################################################################
import sys
import os
import dataset
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
	scores = {}
	final_list = []
	for u in unknowns:
		uRef = u['refNum']
		name = u['name']
		query = "select * from transformations where obs_from=\'%s\' or obs_to=\'%s\'" % (uRef, uRef)
		K = [x for x in db.query(query)]
		S = 0
		if len(K) == 0:
			continue
		for k in K:
			S += k['dmz_err']
		S /= len(K)
		S = 1 / S
		scores[name] = S
		final_list.append(u)
	final_list.sort(key = lambda u : scores[u['name']], reverse = True)

	for u in final_list:
		print(u['name'], u['refNum'], scores[u['name']])

	return final_list




