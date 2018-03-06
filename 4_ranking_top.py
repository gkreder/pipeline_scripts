################################################################################
# Gabe Reder
# gkreder@gmail.com
# Github version
################################################################################
import argparse
import sys
import os
from aux.ranking import rank
################################################################################
# If called directory from the command line
if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--in_db', '-i', dest='in_db', required = True)
	args = parser.parse_args()
	ranked_list = rank(args)
################################################################################