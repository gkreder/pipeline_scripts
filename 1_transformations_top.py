################################################################################
# Gabe Reder
# gkreder@gmail.com
# Github version
################################################################################
import argparse
import sys
import os
from xlrd import open_workbook
from aux.transformations import find_transformations, write_transformations
from aux.merge import merge
import time
################################################################################
def mergedTransformations(args):
	
	args, observations, transformations, adducts, isoforms = find_transformations(args)
	write_transformations(args, observations, transformations, adducts, isoforms)
	merged_data = merge(args, observations, transformations, adducts, isoforms)
	merged_data.writeLog(os.path.join(args.out_dir, '') + 'merged_data.log', args)
	# merged_data.save(os.path.join(args.out_dir, '') + 'merged_data.pkl')
	with open(os.path.join(args.out_dir, '') + 'merged_data.txt', 'w') as f:
		print(merged_data, file = f)
	if os.path.exists(os.path.join(args.out_dir, '') + 'merged_data.sqlite'):
		print('output database found...overwriting')
		os.system('rm %s' % (os.path.join(args.out_dir, '') + 'merged_data.sqlite'))
	merged_data.toSql(os.path.join(args.out_dir, '') + 'merged_data.sqlite')
	print('Done')
	return merged_data
################################################################################
# If called directory from the command line
if __name__ == '__main__':
	trans_tol_mz = 0.01
	adduct_tol_mz = 0.005
	adduct_tol_rt = 0.03
	isoform_tol_mz = 0.0005

	parser = argparse.ArgumentParser()
	parser.add_argument('--trans_file', '-t', required = True, dest = 'trans_file', help = 'Transformation file')
	parser.add_argument('-o', '--out_dir', required = True, dest = 'out_dir', help = 'Output directory')
	parser.add_argument('--trans_tol_mz', default = trans_tol_mz, type = float, help = 'Transformation Tolerance (m/z units) - default %s' % trans_tol_mz)
	parser.add_argument('--adduct_tol_mz', default = adduct_tol_mz, type = float, help = 'Adduct Tolerance (m/z units) - default %s' % adduct_tol_mz)
	parser.add_argument('--adduct_tol_rt', default = adduct_tol_rt, type = float, help = 'Adduct Tolerance (rt units) - default %s' % adduct_tol_rt)
	# parser.add_argument('-i', required = True, dest = 'in_file', help = 'Input File')
	parser.add_argument('--test', dest='test', action='store_true')
	parser.add_argument('--mode', required =True, dest = 'mode', choices = ['pos', 'neg'])
	# parser.add_argument('--input_type', type = str, dest = 'input_type', default = 'meyer')
	parser.add_argument('--no_headers', dest = 'has_header', action='store_false')
	parser.add_argument('--isoform_tol_mz', type = float, default = isoform_tol_mz)
	# parser.add_argument('--knowns_file', '-kf', type = str, required = True, dest = 'knowns_file')
	parser.add_argument('--knowns_db', '-k', type = str, required = True, dest = 'knowns_db')
	parser.add_argument('--input_table_name', '-i', type = str, required = True, dest = 'input_table_name')
	parser.set_defaults(test=False, has_header=True)
	args = parser.parse_args()
	print(args)
	if not os.path.exists(args.out_dir):
		os.system('mkdir ' + args.out_dir)
	args.start_time = time.time()

	mergedTransformations(args)
################################################################################