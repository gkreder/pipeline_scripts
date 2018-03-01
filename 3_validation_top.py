################################################################################
# Gabe Reder - gkreder@gmail.com
# Validation script for checking validity of known-known interactions
# Github version
################################################################################
import sys
import os
import argparse
import dataset
import pickle as pkl
from aux.validation import validate
################################################################################
if __name__ == '__main__':
	# ARGUMENT PARSING
	parser = argparse.ArgumentParser()
	parser.add_argument('--merged_db', '-d', dest='db', type = str, required = True)
	args = parser.parse_args()
	results = validate(args)
