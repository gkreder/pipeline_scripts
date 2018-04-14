################################################################################
import sys
import os
import argparse
################################################################################
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_dir', required = True)
args = parser.parse_args()
################################################################################
dirs = [os.path.join(os.path.abspath(args.input_dir), x) for x in os.listdir(args.input_dir) if '.d' in x and os.path.isdir(os.path.join(os.path.abspath(args.input_dir), x))]

for d in dirs:
	print('msconvert %s' % d)