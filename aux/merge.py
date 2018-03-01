################################################################################
# Pipeline v3.0
# Gabe Reder - gkreder@gmail.com
# 01/2018
################################################################################
import sys
import os
import time
from aux.lib import Obs, MergedSet
################################################################################
def merge(args, observations, transformations, adducts, isoforms):
	start_time = time.time()
	

	# Get known metabolites
	# with open(args.knowns_file, 'r') as f:
	# 	knowns = [x.strip() for x in f]
	# # convert to list of Obs objects
	# knowns = [Obs(x, type = args.input_type) for x in knowns]

	print('\nMerging Data...')
	merged_data = MergedSet.fromObservations(observations)
	merged_data.addAdducts(adducts)
	merged_data.addIsoforms(isoforms)
	merged_data.addTransformations(transformations)
	# knowns = Obs.fromFile(args.knowns_file, args.input_type)
	knowns = Obs.fromTable(args.knowns_db, args.input_table_name, knowns_only = True)
	merged_data.addKnowns(knowns)
	return merged_data