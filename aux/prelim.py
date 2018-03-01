################################################################################
import sys
import os
import dataset
import argparse
import json
from fuzzywuzzy import fuzz, process
from tqdm import tqdm
import time
################################################################################
def scrape_json(json_file):
	scraped_info = {}
	scraped_info['ALIASES'] = []
	scraped_info['MOLECULAR_FORMULA_ALIASES'] = []
	json_data = open(json_file, 'r').read().strip()
	obj = json.loads(json_data)
	obj = obj['Record']
	scraped_info['PUBCHEM_ID'] = obj["RecordNumber"]
	obj = obj['Section']
	for section in obj:
		section_header = section['TOCHeading']
		########################################################################
		# Names and Identifiers Section
		########################################################################
		if section_header == 'Names and Identifiers':
			for subsection in section['Section']:
				subsection_header = subsection['TOCHeading']
				################################################################
				# pubchem name
				################################################################
				if subsection_header == 'Record Title':
					scraped_info['PUBCHEM_NAME'] = subsection['Information'][0]['StringValue']
				elif subsection_header == 'Computed Descriptors':
					name_dicts = [x for x in subsection['Section'] if 'Name' in x['TOCHeading']]
					for n in name_dicts:
						########################################################
						# IUPAC name
						########################################################
						if n['Information'][0]['Name'] == 'IUPAC Name':
							scraped_info['IUPAC_NAME'] = n['Information'][0]['StringValue']
						else:
							scraped_info['ALIASES'].append(n['Information'][0]['StringValue'])
					desc_other = [x for x in subsection['Section'] if 'Name' not in x['TOCHeading']]
					for n in desc_other:
						########################################################
						# InChl
						########################################################
						if n['Information'][0]['Name'] == 'InChI':
							scraped_info['INCHI'] = n['Information'][0]['StringValue']
						########################################################
						# InChl Key
						########################################################
						elif n['Information'][0]['Name'] == 'InChI Key':
							scraped_info['INCHI_KEY'] = n['Information'][0]['StringValue']
						########################################################
						# Canonical SMILES
						########################################################
						elif n['Information'][0]['Name'] == 'Canonical SMILES':
							scraped_info['SMILES'] = n['Information'][0]['StringValue']
						########################################################
						# Any other name/identifier alias
						########################################################
						else:
							scraped_info['ALIASES'].append(n['Information'][0]['StringValue'])
				################################################################
				# molecular formula
				################################################################
				elif subsection_header == 'Molecular Formula':
					for i, formula_entry in enumerate(subsection['Information']):
						if i == 0:
							scraped_info['MOLECULAR_FORMULA'] = formula_entry['StringValue']
						else:
							scraped_info['MOLECULAR_FORMULA_ALIASES'].append(formula_entry['StringValue'])
				################################################################
				# Status (e.g. live vs non-live page)
				################################################################
				elif subsection_header == 'Status':
					for i, status_entry in enumerate(subsection['Information']):
						scraped_info['NON-LIVE'] = status_entry['StringValue']
						print(status_entry['StringValue'])
						sys.exit()

				################################################################
				# CAS Number / other identifiers
				################################################################
				# Just saving CAS number for now, can look here for other 
				# identifers if need - e.g. EC number, etc
				elif subsection_header == 'Other Identifiers':
					identifier_sections = [x for x in subsection['Section']]
					for s in identifier_sections:
						# Find consensus CAS number (one most frequently cited)
						if s['TOCHeading'] == 'CAS':
							found_CAS_numbers = {}
							for entry in s['Information']:
								if entry['StringValue'] in found_CAS_numbers:
									found_CAS_numbers[entry['StringValue']] += 1
								else:
									found_CAS_numbers[entry['StringValue']] = 1
							# The most found CAS number...
							CAS_number = max(found_CAS_numbers.keys(), key=(lambda key: found_CAS_numbers[key]))
							scraped_info['CAS_NUMBER'] = CAS_number
							found_CAS_numbers = list(found_CAS_numbers.keys())
							# found_CAS_numbers.remove(CAS_number)
							scraped_info['CAS_NUMBER_ALIASES'] = found_CAS_numbers
				################################################################
				# Synonyms
				################################################################
				elif subsection_header == 'Synonyms':
					synonym_sections = [x for x in subsection['Section']]
					for syn_sect in synonym_sections:
						# MeSH Entry Terms are exact matches
						if syn_sect['TOCHeading'] == 'MeSH Entry Terms':
							for syn_entry in syn_sect['Information']:
								for s in syn_entry['StringValueList']:
									scraped_info['ALIASES'].append(s)
						########################################################
						# Depositer-supplied synonyms (vs. aliases from pubchem)
						########################################################
						elif syn_sect['TOCHeading'] == 'Depositor-Supplied Synonyms':
							for syn_entry in syn_sect['Information']:
								for s in syn_entry['StringValueList']:
									if 'SYNONYMS' in scraped_info:
										scraped_info['SYNONYMS'].append(s)
									else:
										scraped_info['SYNONYMS'] = [s]
						
		########################################################################
		# Chemical and Physical Properties
		########################################################################
		elif section_header == 'Chemical and Physical Properties':
			subsections = [x for x in section['Section']]
			for subsection in subsections:
				if subsection['TOCHeading'] == 'Computed Properties':
					# Note - just taking the first table here. I only see one on
					# the urea page, just want to be safe. I'm assuming that
					# the first one will always be the pubchem (most validated) one

					# Had a different formatting for one of the entries
					if 'Information' in subsection:
						top_prop_table = subsection['Information'][0]
					else:
						top_prop_table = subsection['Section'][0]
					for row in top_prop_table['Table']['Row']:
						prop_name = [x for x in row['Cell'] if 'StringValue' in x][0]['StringValue']
						########################################################
						# Molecular Weight
						########################################################
						if prop_name == 'Molecular Weight':
							molecular_weight = [x for x in row['Cell'] if 'NumValue' in x][0]['NumValue']
							scraped_info['MOLECULAR_WEIGHT'] = molecular_weight
						########################################################
						# Hydrogen bond donor count
						########################################################
						elif prop_name == 'Hydrogen Bond Donor Count':
							h_bond_donor = [x for x in row['Cell'] if 'NumValue' in x][0]['NumValue']
							scraped_info['H_BOND_DONOR_COUNT'] = h_bond_donor
						########################################################
						# Hydrogen bond acceptor count
						########################################################
						elif prop_name == 'Hydrogen Bond Acceptor Count':
							h_bond_acceptor = [x for x in row['Cell'] if 'NumValue' in x][0]['NumValue']
							scraped_info['H_BOND_ACCEPTOR_COUNT'] = h_bond_acceptor
						########################################################
						# Rotatable bond count
						########################################################
						elif prop_name == 'Rotatable Bond Count':
							rotatable_bond = [x for x in row['Cell'] if 'NumValue' in x][0]['NumValue']
							scraped_info['ROTATABLE_BOND_COUNT'] = rotatable_bond
						########################################################
						# Topological Polar surface area
						########################################################
						elif prop_name == 'Topological Polar Surface Area':
							polar_area = [x for x in row['Cell'] if 'NumValue' in x][0]['NumValue']
							scraped_info['TOPOLOGICAL_POLAR_SURFACE_AREA'] = polar_area
						########################################################
						# xlogp3_aa
						########################################################
						elif prop_name == 'XLogP3-AA':
							xlogp3_aa = [x for x in row['Cell'] if 'NumValue' in x][0]['NumValue']
							scraped_info['XLOGP3_AA'] = xlogp3_aa
							if 'LOGP_CALC' in scraped_info:
								sys.exit('Error: Found XLogP3-AA after already finding XLogP3')
							scraped_info['LOGP_CALC'] = xlogp3_aa
						########################################################
						# xlogp3
						########################################################
						elif prop_name == 'XLogP3':
							xlogp3 = [x for x in row['Cell'] if 'NumValue' in x][0]['NumValue']
							scraped_info['XLOGP3'] = xlogp3
							if 'LOGP_CALC' in scraped_info:
								sys.exit('Error: Found XLogP3 after already finding XLogP3-AA')
							scraped_info['LOGP_CALC'] = xlogp3
						########################################################
						# Formal charge
						########################################################
						elif prop_name == 'Formal Charge':
							formal_charge = [x for x in row['Cell'] if 'NumValue' in x][0]['NumValue']
							scraped_info['FORMAL_CHARGE'] = formal_charge
				

				elif subsection['TOCHeading'] == 'Experimental Properties':
					exp_section = subsection['Section']
					for s in exp_section:
						########################################################
						# Water Solubility
						########################################################
						if s['TOCHeading'] == 'Solubility':
							sol_info = s['Information']
							for sol_entry in sol_info:
								if 'water' in sol_entry['Name'].lower():
									if 'NumValue' in sol_entry.keys():
										if 'WATER_SOLUBILITY' in scraped_info:
											scraped_info['WATER_SOLUBILITY'] += '\t' + str(sol_entry['NumValue']) + ' ' + sol_entry['ValueUnit']
										else:
											scraped_info['WATER_SOLUBILITY'] = str(sol_entry['NumValue']) + ' ' + sol_entry['ValueUnit']
						########################################################
						# Experimental LogP
						########################################################
						elif s['TOCHeading'] == 'LogP':
							logp_info = s['Information']
							for l in logp_info:
								if 'NumValue' in l:
									logp = l['NumValue']
								else:
									logp = l['StringValue']
								if 'LOGP_EXPERIMENTAL' in scraped_info:
									scraped_info['LOGP_EXPERIMENTAL'] += '\t' + str(logp)
								else:
									scraped_info['LOGP_EXPERIMENTAL'] = str(logp)
						########################################################
						# Experimental LogS
						########################################################
						elif s['TOCHeading'] == 'LogS':
							logs_info = s['Information']
							for l in logs_info:
								if 'NumValue' in l:
									logs = l['NumValue']
								else:
									logs = l['StringValue']
								if 'LOGS_EXPERIMENTAL' in scraped_info:
									scraped_info['LOGS_EXPERIMENTAL'] += '\t' + str(logs)
								else:
									scraped_info['LOGS_EXPERIMENTAL'] = str(logs)
						########################################################
						# pKa
						########################################################
						elif s['TOCHeading'] == 'pKa':
							pka_info = s['Information']
							for p in pka_info:
								if 'NumValue' in p:
									pka = str(p['NumValue'])
									if 'ValueUnit' in p:
										pka += ' ' + p['ValueUnit']
								elif 'StringValue' in p:
									pka = p['StringValue']
								else:
									continue
								if 'PKA' in scraped_info:
									scraped_info['PKA'] += '\t' + str(pka)
								else:
									scraped_info['PKA'] = str(pka)
		########################################################################
		# These are unused sections in pubchem that are available
		########################################################################
		elif section_header == '2D Structure':
			continue
		elif section_header == '3D Conformer':
			continue
		elif section_header == 'LCSS':
			continue
		elif section_header == 'Related Records':
			continue
		elif section_header == 'Chemical Vendors':
			continue
		elif section_header == 'Drug and Medication Information':
			continue
		elif section_header == 'Food Additives and Ingredients':
			continue
		elif section_header == 'Agrochemical Information':
			continue
		elif section_header == 'Pharmacology and Biochemistry':
			continue
		elif section_header == 'Use and Manufacturing':
			continue
		elif section_header == 'Identification':
			continue
		elif section_header == 'Safety and Hazards':
			continue
		elif section_header == 'Toxicity':
			continue
		elif section_header == 'Literature':
			continue
		elif section_header == 'Patents':
			continue
		elif section_header == 'Biomolecular Interactions and Pathways':
			continue
		elif section_header == 'Biological Test Results':
			continue
		elif section_header == 'Classification':
			continue
		elif section_header == 'Biologic Description':
			continue
		elif section_header == '3D Status':
			continue
		else:
			print('THIS SHOULDNT PRINT - SECTION UNACCOUNTED FOR')
			print(section_header)
			sys.exit()
		########################################################################
		########################################################################
	return(scraped_info)

def download_pid(pid, args):
	os.system('wget -O %s https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/%s/JSON/' % (os.path.join(args.pubchem_dir,str(pid)+'.json'),str(pid)))
	time.sleep(0.5)


def known_rows(in_file):
	with open(in_file, 'r') as f:
		lines = [x.replace('\n','') for x in f]
	header = lines[0].split('\t')
	lines = lines[1 : ]
	lines = [x.split('\t') for x in lines if x != '']
	lines = [{header[i] : str(l) for i, l in enumerate(x)} for x in lines]
	return lines, header

def yes_or_no(question):
	while "the answer is invalid":
		reply = str(input(question)).lower().strip()
		if reply[:1] == 'o':
			return 'overwrite'
		if reply[:1] == 'e':
			return 'exit'
		if reply[:1] == 'a':
			return 'append'

def add_data(args):
	no_matches = []
	# pubchem_dir = os.path.join('', args.pubchem_dir)
	# json_files = [pubchem_dir + x for x in os.listdir(pubchem_dir) if '.json' in x.lower()]
	if not os.path.exists(os.path.dirname(args.db)):
		os.makedirs(os.path.dirname(args.db))
	db = dataset.connect('sqlite:///' + args.db)
	rows, header = known_rows(args.in_file)
	if 'pid' not in rows[0]:
		sys.exit('Error: The input knowns file must have a column called \"pid\"')
	drop_table = None	
	if args.table_name in db.tables:
		drop_table = yes_or_no('\nTable \"%s\" exists in db what would you like to do? (o)verwrite/(e)xit/(a)ppend: ' % args.table_name)
		if drop_table == 'overwrite':
			print('Overwriting...')
			db[args.table_name].drop()
		elif drop_table == 'exit':
			sys.exit('Exiting...')
		elif drop_table == 'append':
			print('Appending...')
	print('ADDING ROWS')
	formatted_rows = []
	for r in tqdm(rows):
		if drop_table == 'append':
			same_row = [x for x in db.query('select * from %s where compound = \"%s\"' % (args.table_name, r['compound'].replace('\"','')))]
			if len(same_row) > 0:
				db[args.table_name].delete(compound=r['compound'].replace('\"',''))
		pid = r['pid']
		if pid != '':
			# print(pid)
			pid_results = [x for x in db.query('select * from compounds where PUBCHEM_ID = ' + str(pid))]
			if len(pid_results) > 1:
				sys.exit('Error: Multiple PID hits in database for PID = ' + str(pid))
			elif len(pid_results) < 1:
				sys.exit('Error: No PID hits in database for PID = ' + str(pid))
			fid = pid_results[0]['FID']
			r['FID'] = fid
		else:
			no_matches.append(r)
		db[args.table_name].insert(r)
		formatted_rows.append(r)
	if len(no_matches) > 0:
		err_file = os.path.join(os.path.dirname(args.in_file), args.table_name + '_NO_MATCHES.tsv')
		print('Warning: %i compounds had no matching PID in the db' % len(no_matches))
		print('\tWriting these compounds to file %s' % err_file)
		with open(err_file, 'w') as f:
			print('\t'.join(header), file = f)
			for nm in no_matches:
				print('\t'.join([nm[x] for x in header]), file = f)
	# print('\n\nCREATING FORMATTED FILE FOR PIPELINE INPUT')
	# formatted_header = [x for x in db[args.table_name].columns if x != 'id']
	# if not os.path.exists(os.path.dirname(args.formatted_out_file)):
	# 	os.makedirs(os.path.dirname(args.formatted_out_file))
	# with open(args.formatted_out_file, 'w') as f:
	# 	print('\t'.join(formatted_header), file = f)
	# 	for r in tqdm(formatted_rows):
	# 		for x in formatted_header:
	# 			if x not in r.keys():
	# 				r[x] = ''
	# 		print('\t'.join([str(r[col]) for col in formatted_header]), file = f)



def check(args):
	rows, header = known_rows(args.in_file)
	db = dataset.connect('sqlite:///' + args.db)
	pubchem_ids = [str(x.lower().split('.json')[0]) for x in os.listdir(args.pubchem_dir) if '.json' in x.lower()]
	# print([x for x in args.pubchem_dir])
	none = []
	one = []
	more = [] 
	add_pids = []
	for kr in rows:
		# if '\"' in kr['compound']:
			# sys.exit('\n\tFatal: Remove \" character from compound %s' % kr['compound'])
		if '%' in kr['compound']:
			sys.exit('\n\tFatal: Remove %% character from compound %s\n' % kr['compound'])
		if kr['pid'] == '':
			none.append(kr)
			continue
		results = [x for x in db.query('select * from compounds where PUBCHEM_ID = %s' % kr['pid'])]
		if len(results) == 0:
			none.append(kr)
			if str(kr['pid']) in pubchem_ids:
				add_pids.append(str(kr['pid']))
		elif len(results) == 1:
			one.append(kr)
		else:
			more.append(kr)
	print('RESULTS:')
	print('\t%i compounds matched exactly to DB' % len(one))
	print('\t%i compounds have no match' % len(none))
	print('Missing PIDs to add to DB: ')
	for n in none:
		if n['pid'] not in ['', None, '\n']:
			print('%s' % str(n['pid']))
	if args.verbose:
		print('\t\t--------------------------------')
		print('\t\tUnmatched compounds:')
		for cn in none:
			print('\t\t\t%s' % cn['compound'])
		print('\t\t--------------------------------')
	if len(more) > 0:
		print('ERROR: The following %i compounds have multiple matches:' % len(more))
		for m in more:
			print('\t%i, %s' % (m['pid'], m['compound']))
	if len(add_pids) > 0:
		print('ERROR: The following %i PIDS are in the pubchem directory but not in db:' % len(add_pids))
		for p in add_pids:
			print('\t%s' % p)

def add_pids(args):
	db = dataset.connect('sqlite:///%s' % args.db)
	if args.mode == 'file':
		with open(args.pids_list, 'r') as f:
			pids = [x.strip() for x in f]
	elif args.mode == 'string':
		pids = args.pids_list.split(',')
	else:
		sys.exit('Error: add_pids mode must be either file or string')
	print('Downloading PIDS')
	for p in pids:
		download_pid(p, args)
	print('Adding to DB')
	for p in pids:
		jf = os.path.join(args.pubchem_dir, str(p)+'.json')
		jd = scrape_json(jf)
		jd_row = {}
		for x in jd.keys():
			if type(jd[x]) == str:
				jd_row[x] = jd[x]
			elif type(jd[x]) == list:
				jd_row[x] = ';'.join(jd[x])
			elif type(jd[x]) == float:
				jd_row[x] = str(jd[x])
			elif type(jd[x]) == int:
				jd_row[x] = str(jd[x])
			else:
				print(jd)
				print(x)
				print(type(jd[x]))
				sys.exit('Error: unexpected type in scraped data')
		in_db = [x for x in db.query('select * from compounds where PUBCHEM_ID = %s' % str(p))]
		fid = None
		if len(in_db) > 0:
			fid = in_db[0]['FID']
			db['compounds'].delete(PUBCHEM_ID=p)
		if fid in [None, '', 'NULL']:
			fid_max = [x['max(FID)'] for x in db.query('select max(FID) from compounds')][0]
			fid = fid_max + 1
		jd_row['FID'] = fid
		db['compounds'].insert(jd_row)

def download_pids(args):
	if args.mode == 'file':
		with open(args.pids_list, 'r') as f:
			pids = [x.strip() for x in f]
	elif args.mode == 'string':
		pids = args.pids_list.split(',')
	else:
		sys.exit('Error: add_pids mode must be either file or string')
	print('Downloading PIDS')
	for p in pids:
		download_pid(p, args)

