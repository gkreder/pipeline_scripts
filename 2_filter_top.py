################################################################################
# Gabe Reder
# gkreder@gmail.com
# Github version
################################################################################
import argparse
import sys
import os
from aux.filter import filter_transformations
################################################################################
# If called directory from the command line
if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--merged_db', '-i', dest='in_db', required = True)
	parser.add_argument('--out_db', '-o', dest='out_db', required = True)
	args = parser.parse_args()
	filter_transformations(args)
################################################################################