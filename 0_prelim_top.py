################################################################################
# Gabe Reder
# gkreder@gmail.com
# Github version
################################################################################
import sys
import os
import dataset
import argparse
from aux import prelim
################################################################################



if __name__ != '__main__':
	sys.exit()

parser = argparse.ArgumentParser(prog='command')
subparsers = parser.add_subparsers()
subparsers.required = True
subparsers.dest = 'command'

# parser for 'add_data' command
parser_add_prelim= subparsers.add_parser('add_data', help = 'add input data/observations to db')
parser_add_prelim.add_argument('-d', '--db', type = str, required = True)
parser_add_prelim.add_argument('--in_file', '-i', type = str, required = True)
parser_add_prelim.add_argument('--table_name', '-t', type = str, required = True)
# parser_add_prelim.add_argument('-o', '--formatted_out_file', type = str, required = True)

# parser for 'add_pids' command
parser_add_pids = subparsers.add_parser('add_pids', help = 'add input list of pubchem IDs to db and create a formatted input file')
parser_add_pids.add_argument('-d', '--db', type = str, required = True)
parser_add_pids.add_argument('--pubchem_dir', '-p', type = str, required = True)
parser_add_pids.add_argument('-i', '--pids_list', type = str, required = True)
parser_add_pids.add_argument('-m', '--mode', type = str, required = True, choices = ['string', 'file'])

# parser for 'download_pids' command
parser_add_pids = subparsers.add_parser('download_pids', help = 'download pids to json')
parser_add_pids.add_argument('--pubchem_dir', '-p', type = str, required = True)
parser_add_pids.add_argument('-i', '--pids_list', type = str, required = True)
parser_add_pids.add_argument('-m', '--mode', type = str, required = True, choices = ['string', 'file'])

# parser for 'check' command
parser_check = subparsers.add_parser('check', help = 'preliminary Check on input list of knowns')
parser_check.add_argument('-d', '--db', type = str, required = True)
parser_check.add_argument('--in_file', '-i', type = str, required = True)
parser_check.add_argument('--pubchem_dir', '-p', type = str, required = True)
parser_check.add_argument('-v', '--verbose', action='store_true')

# # parser for 'match_input' command
# parser_check = subparsers.add_parser('match_input', help = 'add ID labels to input file')
# parser_check.add_argument('-d', '--db', type = str, required = True)
# parser_check.add_argument('--knowns_file', '-k', type = str, required = True)
# parser_check.add_argument('-o', '--out_file', required = True)


args = parser.parse_args()

if args.command == 'add_data':
	prelim.add_data(args)
elif args.command == 'add_pids':
	prelim.add_pids(args)
elif args.command == 'check':
	prelim.check(args)
# elif args.command == 'match_input':
# 	prelim.match_input(args)
elif args.command == 'download_pids':
	prelim.download_pids(args)
else:
	sys.exit('Error: Unrecognized command - %s' % args.command)