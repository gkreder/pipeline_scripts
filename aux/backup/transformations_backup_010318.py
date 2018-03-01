################################################################################
# Gabe Reder - gkreder@gmail.com
# 08/2017
# 
# Finds candidate relationships between metabolites 
################################################################################
import sys
import os
import time
from xlrd import open_workbook
from aux.lib import Obs, InTrans, OutTrans, check_time
def find_transformations(args):
	start_time = time.time()
	############################################################################
	# Get lines from input file and filter/save
	##########################################################
	transformations = InTrans.fromExcel(args.trans_file)
	max_tdiff = max([abs(x.dmz) for x in transformations])

	# open input file and read lines (ignoring NA lines)
	with open(args.in_file, 'r') as f:
		# Scrape text lines
		lines = [line.strip() for line in f]
	# Throw away header
	lines = lines[1 : ]
	# Convert lines to Obs objects
	lines = [Obs(x) for x in lines]
	# sort lines
	lines = sorted(lines, key = lambda line : line.mz)
	# Check if all the same method
	mc = set([x.method for x in lines])
	if len(mc) > 1:
		sys.exit('''Error: Multiple methods in the same input file.
			Must separate inputs from different methods''')
	if args.test:
		lines = lines[0 : 500]
	############################################################################
	# Find Transformations
	############################################################################
	saved_transformations = []
	j_start = 0
	# method_count = {}
	for i, line_outer in enumerate(lines):
		check_time(i, lines, start_time)
		# if line_outer.method in method_counts:
			# method_counts[line_outer.method] = method_counts[method] + 1
		# else:
			# method_counts[line_outer.method] = 1

		# Indexing resets itself to avoid unecesarry comparisons
		j_temp = j_start
		for j in range(j_start, len(lines)):
			line_inner = lines[j]
			delta_mz = line_inner.mz - line_outer.mz
			if abs(delta_mz) > abs(max_tdiff) + args.trans_tol:
				if j >= i:
					break
				else:
					if j > j_start:
						j_temp = j
					continue
			if line_outer.method == line_inner.method:
				# temp method
				def cm(x):
					m1 = abs(x.dmz - delta_mz)
					m2 = abs(x.dmz + delta_mz)
					return(min(m1, m2))
				# given our transformations, what's the closest possible one to
				# the one we've observed
				closest_match = min(transformations, key = lambda x : cm(x))
				dist_current = abs(closest_match.dmz - delta_mz)
				dist_reverse = abs(closest_match.dmz + delta_mz)
				md = min(dist_reverse, dist_current)
				# Check if within trans tolerance
				if md <= args.trans_tol:
					obs_to = line_inner
					obs_from = line_outer
					# Make sure the order of transformation is correct
					if dist_reverse < dist_current:
						obs_to = line_outer
						obs_From = line_inner
					t_o = OutTrans(closest_match, obs_from, obs_to)
					saved_transformations.append(t_o)
		j_start = j_temp
	############################################################################
	# Write output file / return output
	############################################################################
	out_dir = os.path.join(args.out_dir, '')
	if not os.path.exists(out_dir + 'transformations/'):
		os.system('mkdir ' + out_dir + 'transformations/')
	t_log = out_dir + 'transformations/transformations.log'
	t_out = out_dir + 'transformations/transformations.tsv'
	with open(t_out, 'w') as f:
		print(OutTrans.header(), file = f)
		for t_o in saved_transformations:
			print(t_o.line(), file = f)
	end_time = time.time()
	with open(t_log, 'w') as f:
		print(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()), file = f)
		print('Runtime - ' + str(end_time - start_time), file = f)
		print('Infile - ' + args.in_file, file = f)
		print('Tolerance - ' + str(args.trans_tol), file = f)
		print('Input observations - ' + str(len(lines)), file = f)
		print('Found transformations - ' + str(len(saved_transformations)), file = f)
	return saved_transformations