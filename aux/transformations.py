################################################################################
# Pipeline v2.0
# Gabe Reder - gkreder@gmail.com
# 01/2018
# 
# Finds candidate relationships between metabolites 
################################################################################
import sys
import os
import time
from xlrd import open_workbook
from aux.lib import Obs, RefTrans, Trans, RefAdduct, RefAdductTrans, AdductTrans, Isoform
from tqdm import tqdm
import dataset
################################################################################
def find_transformations(args):
	############################################################################
	# Get lines from input file and filter/save
	##########################################################
	db = dataset.connect('sqlite:///' + os.path.abspath(args.knowns_db))
	if args.input_table_name not in db.tables:
		sys.exit(print('\nError: No table named %s found in database %s\n' % (args.input_table_name, args.knowns_db)))

	transformations_ref = RefTrans.fromExcel(args.trans_file)
	max_tdiff = max([abs(x.dmz) for x in transformations_ref])

	adductTrans_ref = RefAdductTrans.fromExcel(args.trans_file, args.mode)
	max_adiff = max([abs(x.dmz) for x in adductTrans_ref])


	# # open input file and read lines (ignoring NA lines)
	# with open(args.in_file, 'r') as f:
	# 	# Scrape text lines
	# 	lines = [line.strip() for line in f]
	# # Throw away header
	# lines = lines[1 : ]
	# # Convert lines to Obs objects
	# lines = [Obs(x, input_type = args.input_type) for x in lines]

	# get observations from input file
	# lines = Obs.fromFile(args.in_file, input_type = args.input_type)
	lines = Obs.fromTable(args.knowns_db, args.input_table_name)

	# sort lines
	lines_trans = sorted(lines, key = lambda line : float(line.mz))
	lines_adducts = sorted(lines, key = lambda line : float(line.rt))
	# Check if all the same method
	# mc = set([x.method for x in lines_trans])
	# if len(mc) > 1 and args.input_type != 'shuo_sql':
	# 	sys.exit('''Error: Multiple methods in the same input file.
	# 		Must separate inputs from different methods''')
	if args.test:
		lines_trans = lines_trans[0 : 500]
		lines_adducts = lines_adducts[0 : 500]
		lines = lines[0 : 500]


	# Check that there are no repeats in the input data
	print('Checking input data for repeats...')
	bad_repeats = []
	warning_repeats = []
	for i,o1 in enumerate(tqdm(lines_trans)):
		for j in range(i, len(lines_trans)):
			o2 = lines_trans[j]
			if i == j:
				continue
			if o2.mz > o1.mz:
				break
			if o1.name != o2.name:
				continue
			if o1.rt != o2.rt:
				warning_repeats.append(o1)
			bad_repeats.append(o1.name)
	if len(bad_repeats) > 0:
		for br in bad_repeats:
			print(br)
		sys.exit('ERROR (transformations.py): Found repeat duplicate entries for the compounds above in the input data table - fix this before proceeding')

	if len(warning_repeats) > 0:
		print('---------------------------------------------------------------')
		print('WARNING - The Following compounds have suspect duplicate entries (with different RTs)')
		for wr in warning_repeats:
			print('\t%s' % wr)
		print('---------------------------------------------------------------')
	############################################################################
	# Find Transformations
	############################################################################
	print('\nFinding Transformations...')
	redun_transformations_from = {} # keys = names, values = list of pids
	redun_transformations_to = {}
	saved_transformations = []
	saved_observations = []
	saved_observations_indices = []
	saved_isoforms = []
	j_start = 0
	# method_count = {}
	# for i, line_outer in enumerate(lines_trans):
	for i, line_outer in enumerate(tqdm(lines_trans)):
		if i not in saved_observations_indices:
			saved_observations_indices.append(i)
			saved_observations.append(line_outer)
		# check_time(i, lines_trans, args.start_time, tabs = 1)
		# Indexing resets itself to avoid unecesarry comparisons
		j_temp = j_start
		for j in range(i + 1, len(lines_trans)):
			line_inner = lines_trans[j]
			if j not in saved_observations_indices:
				saved_observations_indices.append(j)
				saved_observations.append(line_inner)
			delta_mz = float(line_inner.mz) - float(line_outer.mz)
			if abs(delta_mz) > abs(max_tdiff) + args.trans_tol_mz:
				if j >= i:
					break
				else:
					if j > j_start:
						j_temp = j
					continue
			if i == j:
				continue

			# if line_outer.method == line_inner.method:
			# temp method
			def cm(x):
				m1 = abs(x.dmz - delta_mz)
				m2 = abs(x.dmz + delta_mz)
				return(min(m1, m2))
			# given our transformations, what's the closest possible one to
			# the one we've observed
			closest_match = min(transformations_ref, key = lambda x : cm(x))
			dist_current = abs(closest_match.dmz - delta_mz)
			dist_reverse = abs(closest_match.dmz + delta_mz)
			md = min(dist_reverse, dist_current)
			# Check if within trans tolerance				
			if closest_match.name.lower() == 'isoform':
				if md <= args.isoform_tol_mz:
					obs_from = line_inner
					obs_to = line_outer
					i_s = Isoform(obs_from, obs_to)
					saved_isoforms.append(i_s)
			else:
				if md <= args.trans_tol_mz:
					obs_to = line_inner
					obs_from = line_outer


					# Make sure the order of transformation is correct
					# This was a bug, but why did it make the accuracy better?
					if dist_reverse < dist_current:
						obs_to = line_outer
						obs_from = line_inner
					t_s = Trans(closest_match, obs_from, obs_to)
					if not args.no_filter:
						# Don't allow redundant PID transformations
						if obs_from.known:
							if obs_to.name in redun_transformations_to:
								if obs_from.pid in redun_transformations_to[obs_to.name]:
									continue
								else:
									redun_transformations_to[obs_to.name].append(obs_from.pid)
							else:
								redun_transformations_to[obs_to.name] = [obs_from.pid]
						if obs_to.known:
							if obs_from.name in redun_transformations_from:
								if obs_to.pid in redun_transformations_from[obs_from.name]:
									continue
								else:
									redun_transformations_from[obs_from.name].append(obs_to.pid)
							else:
								redun_transformations_from[obs_from.name] = [obs_to.pid]

						# Filter impossible transformations
						


					if obs_from == obs_to:
						print(line_outer)
						print(line_inner)
						sys.exit('Error (transformations.py) - the two observations should not be the same')
					saved_transformations.append(t_s)
		j_start = j_temp

	# Check for duplicates
	# for i_outer, t_outer in enumerate((saved_transformations)):
	# 	for j_outer in range(i_outer + 1, len(saved_transformations)):
	# 		t_inner = saved_transformations[j_outer]
	# 		if t_outer.name == t_inner.name and t_outer.trans.name == t_inner.trans.name:
	# 			print(t_inner, t_outer)
	############################################################################
	# Find Adducts - see commens in find transformations section
	############################################################################
	print('\nFinding Adducts...')
	saved_adducts = []
	j_start = 0
	# for i, line_outer in enumerate(lines_adducts):
	for i, line_outer in enumerate(tqdm(lines_adducts)):
		# check_time(i, lines_adducts, args.start_time, tabs = 1)
		j_temp = j_start
		for j in range(j_start, len(lines_adducts)):
			line_inner = lines_adducts[j]
			delta_rt = float(line_inner.rt) - float(line_outer.rt)
			if abs(delta_rt) > args.adduct_tol_rt:
				if j >= i:
					break
				else:
					if j > j_start:
						j_temp = j
					continue
			if i == j:
				continue
			if line_outer.method == line_inner.method:
				delta_mz = float(line_inner.mz) - float(line_outer.mz)
				def cm(x):
					m1 = abs(x.dmz - delta_mz)
					m2 = abs(x.dmz + delta_mz)
					return(min(m1, m2))
				closest_match = min(adductTrans_ref, key = lambda x : cm(x))
				dist_current = abs(closest_match.dmz - delta_mz)
				dist_reverse = abs(closest_match.dmz - delta_mz)
				md = min(dist_current, dist_reverse)
				if md <= args.adduct_tol_mz:
					line_to = line_inner
					line_from = line_outer
					if dist_reverse < dist_current:
						line_to = line_outer
						line_from = line_inner
					a_s = AdductTrans(closest_match, line_from, line_to)
					saved_adducts.append(a_s)
		j_start = j_temp
	args.t_end_time = time.time()
	return args, saved_observations, saved_transformations, saved_adducts, saved_isoforms

def write_transformations(args, saved_observations, saved_transformations, saved_adducts, saved_isoforms):
	############################################################################
	# Write output file / return output
	############################################################################
	out_dir = os.path.join(args.out_dir, '')
	if not os.path.exists(out_dir + 'merged_intermediate_files/'):
		os.system('mkdir ' + out_dir + 'intermediate_files/')
	t_log = out_dir + 'intermediate_files/transformations.log'
	t_out = out_dir + 'intermediate_files/transformations.tsv'
	a_log = out_dir + 'intermediate_files/adducts.log'
	a_out = out_dir + 'intermediate_files/adducts.tsv'
	i_log = out_dir + 'intermediate_files/isoforms.log'
	i_out = out_dir + 'intermediate_files/isoforms.tsv'
	end_time = time.time()
	runtime = str(end_time - args.start_time)

	Trans.write_out(t_out, saved_transformations)
	Trans.write_log(t_log, runtime, args, saved_observations,saved_transformations)

	AdductTrans.write_out(a_out, saved_adducts)
	AdductTrans.write_log(a_log, runtime, args, saved_observations, saved_adducts)

	Isoform.write_out(i_out, saved_isoforms)
	Isoform.write_log(i_log, runtime, args, saved_observations, saved_isoforms)

	