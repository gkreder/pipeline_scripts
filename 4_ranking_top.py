################################################################################
# Gabe Reder
# gkreder@gmail.com
# Github version
################################################################################
import argparse
import sys
import os
from aux.ranking import clean_rank
################################################################################
# If called directory from the command line
if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--in_db', '-i', dest='in_db', required = True)
	args = parser.parse_args()
	final_list = clean_rank(args)

	# return ranked_list
################################################################################