################################################################################
# Gabe Reder - gkreder@gmail.com
# 01/2018
################################################################################
import pickle as pkl
import os
import math
import time
import sys
from xlrd import open_workbook
import copy
import dataset
import numpy as np
from tqdm import tqdm
from rdkit import Chem
import yaml
################################################################################
# to sql row
################################################################################
def toSqlRow(o):
	n_row = {}
	for var, value in vars(o).items():
		if value == None:
			n_row[var] = None
		elif type(value) == str:
			n_row[var] = value
		elif type(value) in [float, int, bool]:
			n_row[var] = value
		elif type(value) in [dict]:
			n_row[var] = str(value)
		elif type(value) == DSet:
			n_row[var] = ';'.join(str(x.refNum) for x in value.data)
		elif type(value) == list or type(value) == set:
			n_row[var] = ';'.join([str(x.refNum) for x in value])
		elif type(value) in [RefAdduct]:
			n_row[var] = value.name
		elif type(value) in [RefTrans, RefAdductTrans]:
			n_row[var] = value.name
		elif type(value) in [Node, Obs, Edge, Trans]:
			n_row[var] = value.refNum
		else:
			print(vars(o).items(), var, value)
			print('')
			print(type(value))
			sys.exit('\nError: Unrecognized object attribute type in sql row conversion')
	return n_row

################################################################################
# Objects/Classes
################################################################################
class Obs:
	"Class to define observations from input files"
	
	# Assuming from input template file
	def __init__(self, knowns_db, table_name, index):
		db = dataset.connect('sqlite:///' + knowns_db)
		db_results = [x for x in db.query('select * from %s where id = %s' % (table_name, str(index)))]
		if len(db_results) > 1:
			sys.exit('Error: Found multiple index hits for index %s' % index)
		if len(db_results) == 0:
			sys.exit('Error: No PID hits for index %s' % index)
		db_hit = db_results[0]
		if 'mz' in db_hit:
			self.mz = db_hit['mz']
		elif 'mz_calc' in db_hit:
			self.mz = db_hit['mz_calc']
		elif 'mz_obs' in db_hit:
			self.mz = db_hit['mz_obs']
		else:
			sys.exit('Error: No mz entry found in db line for index %s' % index)

		if 'rt' in db_hit:
			self.rt = db_hit['rt']
		elif 'ri' in db_hit:
			self.rt = db_hit['ri']
		else:
			sys.exit('Error: No rt entry found in db line for index %s' % index)

		self.method = db_hit['method']
		self.name = db_hit['compound']
		self.refNum = None
		if 'pid' in db_hit and db_hit['pid'] not in ['None', None, 'NULL', '']:
			if 'FID' not in db_hit or db_hit['FID'] in ['None', None, 'NULL', '']:
				sys.exit('\nError: index %s in table name %s has pid entry but no FID entry' % (str(index), table_name))
			self.pid = db_hit['pid']
			self.FID = db_hit['FID']
			self.known = True
		else:
			self.known = False
			self.pid = None
			self.FID = None

		if self.known:
			compound_hits = [x for x in db.query('select * from compounds where FID = %s' % str(self.FID))]
			if len(compound_hits) == 0:
				sys.exit('\nError: no match in compounds table for FID %s' % str(self.FID))
			elif len(compound_hits) > 1:
				sys.exit('\nError: multiple matches in compounds table for FID %s' % str(self.FID))
			c_hit = compound_hits[0]
			for var, val in c_hit.items():
				setattr(self, var, val)

	def __str__(self):
		return self.name

	# takes in file and returns list of observations in file
	def fromFile(fname, input_type, header = True):
		with open(fname, 'r') as f:
			lines = [x.strip() for x in f]
		if header:
			lines = lines[1 : ]

		obs_out = [Obs(x, input_type=input_type) for x in lines]
		return obs_out

	def fromTable(db_name, table_name, knowns_only = False):
		db = dataset.connect('sqlite:///' + db_name)
		if knowns_only:
			indices = [x['id'] for x in db.query('select distinct id from %s where pid != \'\' and not(pid isnull)' % table_name)]
		else:
			indices = [x['id'] for x in db.query('select distinct id from %s' % table_name)]
		observations = [Obs(db_name, table_name, i) for i in indices]
		return observations

	def fromSqlRow(row):
		# row = [x for x in db.query('select * from observations where refNum = ' + str(refNum))][0]
		o = Obs([], input_type = 'sql')
		for var, val in row.items():
			setattr(o, var, val)
		return o

	def equals(self, o):
		e = True
		if self.mz != o.mz:
			e = False
		if self.rt != o.rt:
			e = False
		if self.method != o.method:
			e = False
		# if self.metabolite != o.metabolite:
			# e = False
		# if self.metab != o.metab:
			# e = False
		if self.name != o.name:
			e = False
		return e

	# os is either a list or a set (some iterable) of observations
	def isIn(self, os):
		for o in os:
			if self.equals(o):
				return True
		return False

	def isKnown(self):
		return self.known


class RefTrans:
	"Class to define input (reference) transformations"
	def __init__(self, name, mz, atom_dict):
		self.name = name
		self.dmz = float(mz)
		self.refNum = None
		self.atoms = atom_dict

	# get transformations from excel file
	def fromExcel(e_file):
		has_isoforms = False
		wb = open_workbook(e_file)
		transformation_sheet = wb.sheet_by_name('Common chemical relationships')
		mz_list = [x.value for x in transformation_sheet.col(7)[2 : ]]
		names = [x.value for x in transformation_sheet.col(0)[2 : ]]
		delta_H_vals = [x.value for x in transformation_sheet.col(1)[2 : ]]
		delta_C_vals = [x.value for x in transformation_sheet.col(2)[2 : ]]
		delta_N_vals = [x.value for x in transformation_sheet.col(3)[2 : ]]
		delta_O_vals = [x.value for x in transformation_sheet.col(4)[2 : ]]
		delta_P_vals = [x.value for x in transformation_sheet.col(5)[2 : ]]
		delta_S_vals = [x.value for x in transformation_sheet.col(6)[2 : ]]

		for n in names:
			if n.lower() == 'isoform':
				has_isoforms = True
		transformations = list(zip(mz_list, names))
		for i, (m,n) in enumerate(transformations):
			atom_dict = {
			'H' : int(delta_H_vals[i]),
			'C' : int(delta_C_vals[i]),
			'N' : int(delta_N_vals[i]),
			'O' : int(delta_O_vals[i]),
			'P' : int(delta_P_vals[i]),
			'S' : int(delta_S_vals[i]),
				}
			transformations[i] = (m,n,atom_dict)

		transformations = [RefTrans(n, m, a) for (m,n,a) in transformations]
		if not has_isoforms:
			sys.exit('Error: please add isoform to chemical trasnformations sheet')
		return transformations

	def fromSqlRow(row):
		name = row['name']
		dmz = row['dmz']
		rt = RefTrans(name, dmz)
		rt.refNum = row['refNum']
		return rt


class Trans:
	"Class to define output (proposed) transformations"
	def __init__(self, closest_match, obs_from, obs_to):
		delta_mz_obs = float(obs_to.mz) - float(obs_from.mz)
		# proposed transformation (Trans object)
		self.trans = closest_match
		# delta_mz observed
		self.dmz_obs = delta_mz_obs
		# delta_mz error
		self.dmz_err = abs(delta_mz_obs - closest_match.dmz)
		# metabolite from
		self.obs_from = obs_from
		# metabolite to
		self.obs_to = obs_to
		self.refNum = None
		# method
		self.method = obs_from.method
		self.name = obs_from.name + '--->' + obs_to.name
		if obs_to.method != self.method:
			sys.exit('Error: Different observation methods for proposed transformation')
		
		# For now, creating a single attribute (known_correct) that encapsulates
		# two conditions. For known_correct to be True, must have a trans that is
		# both known--->known and also that the proposed transformation between
		# the two knowns is correct
		# EDIT - I'm also adding in another attribute called known_known
		# because I might want to distinguish between things that we are 
		# sure are incorrect (that must be known--->known)
		if self.obs_from.known and self.obs_to.known:
			self.known_known = True
			if Trans.correct(closest_match, obs_from, obs_to):
				self.known_correct = True
			else:
				self.known_correct = False
		else:
			self.known_known = False
			self.known_correct = False

	def __str__(self):
		s = 'Transformation: ' + str(self.obs_from) + '--->' + str(self.obs_to)
		return s

	def correct(closest_match, obs_from, obs_to):
		form_from = obs_from.SMILES
		form_to = obs_to.SMILES
		v_from = atomVector(form_from)
		v_to = atomVector(form_to)
		v_diff = vecDifference(v_from, v_to)
		v_expected = closest_match.atoms
		atoms_observed = v_diff
		atoms_expected = v_expected
		if vecEquals(v_diff, v_expected):
			return True
		else:
			return False

	# Return tsv line given a Trans object
	def line(self, delimiter = '\t'):
		l = []
		l.append(self.obs_from.name + '--->' + self.obs_to.name)
		l.append(self.trans.name)
		l.append(str(self.dmz_obs))
		l.append(str(self.trans.dmz - self.dmz_obs))
		l.append(self.obs_from.name)
		l.append(str(self.obs_from.mz))
		l.append(str(self.obs_from.rt))
		l.append(self.obs_to.name)
		l.append(str(self.obs_to.mz))
		l.append(str(self.obs_to.rt))
		l.append(self.method)
		return(delimiter.join(l))
	# Return a tsv header
	def header(delimiter = '\t'):
		h = []
		h.append('transformation')
		h.append('closest_reaction')
		h.append('delta_mz_observed')
		h.append('delta_mz_error - (actual - observed)')
		h.append('metabolite_from')
		h.append('mz_from')
		h.append('rt_from')
		h.append('metabolite_to')
		h.append('mz_to')
		h.append('rt_to')
		h.append('method')
		return(delimiter.join(h))

	def write_log(t_log, runtime, args, lines, saved_transformations):
		with open(t_log, 'w') as f:
			print(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()), file = f)
			print('Runtime - ' + runtime, file = f)
			print('Input DB / Table - %s, %s' % (args.knowns_db, args.input_table_name), file = f)
			print('Tolerance - ' + str(args.trans_tol_mz), file = f)
			print('Input observations - ' + str(len(lines)), file = f)
			print('Found transformations - ' + str(len(saved_transformations)), file = f)

	def write_out(t_out, saved_transformations):
		if not os.path.exists(t_out):
			os.system('touch ' + t_out)
		with open(t_out, 'w') as f:
			print(Trans.header(), file = f)
			for t_s in saved_transformations:
				print(t_s.line(), file = f)

	def fromSqlRow(row, observations_dict, refTransformations_dict):
		obs_from = observations_dict[row['obs_from_refNum']]
		obs_to = observations_dict[row['obs_to_refNum']]
		closest_match = refTransformations_dict[row['refTransformation_refNum']]
		trans = Trans(closest_match, obs_from, obs_to)
		trans.refNum = row['refNum']
		return trans

	def known_known(self):
		if obs_from.isKnown and obs_to.isKnown:
			return True
		else:
			return False

class RefAdduct:
	"Class to define the input adducts"
	def __init__(self, name, dmz, mode):
		self.name = name
		self.dmz = dmz
		self.mode = mode
	def posAdducts(adducts):
		pos = [x for x in adducts if 'pos' in x.mode.lower() or 'all' in x.mode.lower()]
		return pos
	def negAdducts(adducts):
		neg = [x for x in adducts if 'neg' in x.mode.lower() or 'all' in x.mode.lower()]
		return neg
class RefAdductTrans:
	"Transformation BETWEEN two reference adducts"
	def __init__(self, add_from, add_to, mode, m_plus_1=False):
		# Not referring plainly to M--->M+1
		if not m_plus_1:
			self.add_from = add_from
			self.add_to = add_to
		else:
			m_adduct = RefAdduct('M', 0.0, mode)
			m_plus_adduct = RefAdduct('M+1', 1.0033548378, 'all')
			self.add_from = m_adduct
			self.add_to = m_plus_adduct
		self.dmz = self.add_to.dmz - self.add_from.dmz
		self.mode = mode
		self.name = self.add_from.name + '--->' + self.add_to.name
		self.refNum = None
	def fromExcel(e_file, mode):
		wb = open_workbook(e_file)
		transformation_sheet = wb.sheet_by_name('Common adducts')
		# Save mz vals and adduct names (removing blank rows)
		names = [x.value for x in transformation_sheet.col(0)[1 : ] if x.value != '']
		mz_list = [x.value for x in transformation_sheet.col(1)[1 : ] if x.value != '']
		modes = [x.value for x in transformation_sheet.col(2)[1 : ] if x.value != '']
		adducts = list(zip(names, mz_list, modes))
		adducts = [RefAdduct(name, dmz, mode) for (name, dmz, mode) in adducts]
		if mode == 'pos':
			adducts = RefAdduct.posAdducts(adducts)
		elif mode == 'neg':
			adducts = RefAdduct.negAdducts(adducts)
		else:
			sys.exit('Error: Unrecognized mode entered')
		adductTrans_list = []
		for i, add_i in enumerate(adducts):
			if add_i.name.upper() == 'M+1':
				a_t = RefAdductTrans(None, None, mode, m_plus_1=True)
				adductTrans_list.append(a_t)
			for j, add_j in enumerate(adducts):
				if i != j:
					a_t = RefAdductTrans(add_i, add_j, mode)
					adductTrans_list.append(a_t)
		return adductTrans_list

	def fromSqlRow(row):
		a = RefAdduct([], [], '')
		a.dmz = row['dmz']
		a.mode = row['mode']
		a.name = row['name']
		a.refNum = row['refNum']
		return a




class AdductTrans:
	"Proposed found transformation corresponding to two adducts of a metabolite"
	def __init__(self, closest_match, obs_from, obs_to):
		delta_mz_obs = float(obs_to.mz) - float(obs_from.mz)
		self.adductTrans = closest_match
		self.obs_from = obs_from
		self.obs_to = obs_to
		self.dmz_obs = delta_mz_obs
		self.method = obs_from.method
		self.name = obs_from.name + '--->' + obs_to.name
		if obs_to.method != self.method:
			sys.exit('Error: Different observation methods for proposed adduct transformation')
		if type(closest_match) == RefAdductTrans:
			self.mode = closest_match.mode
		else:
			self.mode = self.method
	def __str__(self):
		s = 'Adduct: ' + str(self.obs_from) + '--->' + str(self.obs_to)
		return s
	# Return tsv line given a PropTrans object
	def line(self, delimiter = '\t'):
		l = []
		l.append(self.obs_from.name + '--->' + self.obs_to.name)
		l.append(self.adductTrans.name)
		l.append(str(self.dmz_obs))
		l.append(str(self.adductTrans.dmz - self.dmz_obs))
		l.append(self.obs_from.name)
		l.append(str(self.obs_from.mz))
		l.append(str(self.obs_from.rt))
		l.append(self.obs_to.name)
		l.append(str(self.obs_to.mz))
		l.append(str(self.obs_to.rt))
		l.append(self.method)
		return(delimiter.join(l))
	# Return a tsv header
	def header(delimiter = '\t'):
		h = []
		h.append('transformation')
		h.append('closest_adductTrans')
		h.append('delta_mz_observed')
		h.append('delta_mz_error - (actual - observed)')
		h.append('metabolite_from')
		h.append('mz_from')
		h.append('rt_from')
		h.append('metabolite_to')
		h.append('mz_to')
		h.append('rt_to')
		h.append('method')
		return(delimiter.join(h))
	def write_out(a_out, saved_adducts):
		with open(a_out, 'w') as f:
			print(AdductTrans.header(), file = f)
			for a_s in saved_adducts:
				print(a_s.line(), file = f)
	def write_log(a_log, runtime, args, lines, saved_adducts):
		with open(a_log, 'w') as f:
			print(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()), file = f)
			print('Runtime - ' + runtime, file = f)
			print('Input DB / Table - %s, %s' % (args.knowns_db, args.input_table_name), file = f)
			print('M/Z Tolerance - ' + str(args.adduct_tol_mz), file = f)
			print('RT Tolerance - ' + str(args.adduct_tol_rt), file = f)
			print('Input observations - ' + str(len(lines)), file = f)
			print('Found Adducts - ' + str(len(saved_adducts)), file = f)


	def fromSqlRow(row, observations_dict, refAdducts_dict):
		obs_from = observations_dict[row['obs_from_refNum']]
		obs_to = observations_dict[row['obs_to_refNum']]
		closest_match = refAdducts_dict[row['refAdduct_refNum']]
		adduct = AdductTrans(closest_match, obs_from, obs_to)
		adduct.refNum = row['refNum']
		return adduct


class Isoform():
	"Isoform meaning two observations that look the same"
	def __init__(self, obs_from, obs_to):
		self.obs_from = obs_from
		self.obs_to = obs_to
		self.dmz_obs = float(obs_to.mz) - float(obs_from.mz)
		self.drt = float(obs_to.rt) - float(obs_from.rt)
		self.method = obs_from.method
		self.name = obs_from.name + '--->' + obs_to.name
		if obs_to.method != self.method:
			sys.exit('Error: Different observation methods for proposed isoform')

	def __str__(self):
		s = 'Isoform: ' + str(self.obs_from) + '--->' + str(self.obs_to)
		return s

	# Return tsv line given a Isoform object
	def line(self, delimiter = '\t'):
		l = []
		l.append(self.obs_from.name + '--->' + self.obs_to.name)
		l.append(str(self.dmz_obs))
		l.append(str(0.0 - self.dmz_obs))
		l.append(self.obs_from.name)
		l.append(str(self.obs_from.mz))
		l.append(str(self.obs_from.rt))
		l.append(self.obs_to.name)
		l.append(str(self.obs_to.mz))
		l.append(str(self.obs_to.rt))
		l.append(self.method)
		return(delimiter.join(l))
	# Return a tsv header
	def header(delimiter = '\t'):
		h = []
		h.append('transformation')
		h.append('delta_mz_observed')
		h.append('delta_mz_error - (actual - observed)')
		h.append('metabolite_from')
		h.append('mz_from')
		h.append('rt_from')
		h.append('metabolite_to')
		h.append('mz_to')
		h.append('rt_to')
		h.append('method')
		return(delimiter.join(h))
	def write_out(i_out, saved_isoforms):
		with open(i_out, 'w') as f:
			print(Isoform.header(), file = f)
			for i_s in saved_isoforms:
				print(i_s.line(), file = f)
	def write_log(i_log, runtime, args, lines, saved_isoforms):
		with open(i_log, 'w') as f:
			print(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()), file = f)
			print('Runtime - ' + runtime, file = f)
			print('Input DB / Table - %s, %s' % (args.knowns_db, args.input_table_name), file = f)
			print('M/Z Tolerance - ' + str(args.adduct_tol_mz), file = f)
			print('RT Tolerance - ' + str(args.adduct_tol_rt), file = f)
			print('Input observations - ' + str(len(lines)), file = f)
			print('Found Isoforms - ' + str(len(saved_isoforms)), file = f)

	def fromSqlRow(row, observations_dict):
		# print(observations_dict)
		# print(row['obs_from_refNum'])
		# # print(row['obs_to_refNum'])
		obs_from = observations_dict[row['obs_from_refNum']]
		obs_to = observations_dict[row['obs_to_refNum']]
		i = Isoform(obs_from, obs_to)
		i.refNum = row['refNum']
		return(i)

class DSet:
	def __init__(self):
		self.observations = set([])
		self.data = set([])
	def add(self, d):
		o_from = d.obs_from
		o_to = d.obs_to
		self.observations.add(o_from)
		self.observations.add(o_to)
		self.data.add(d)
	def contains(self, d):
		o_from = d.obs_from
		o_to = d.obs_to
		if o_from in self.observations or o_to in self.observations:
			return True
		else:
			return False
	def union(self, ds):
		ds_out = DSet()
		ds_out.observations = self.observations.union(ds.observations)
		ds_out.data = self.data.union(ds.data)
		return ds_out

	def difference(self, ds):
		ds_out = Dset()
		ds_out.observations = self.observations.difference(ds.observations)
		ds_out.data = self.data.difference(ds.data)
		return ds_out

	def remove(self, d):
		self.data.remove(d)
		self.observations.remove(d.obs_from)
		self.observations.remove(d.obs_to)

	def fromList(d_list):
		ds = DSet()
		for d in d_list:
			ds.add(d)
		return ds



class Edge:
	"connecting two nodes"
	def __init__(self, node_from, node_to, transformation):
		self.node_from = node_from
		self.node_to = node_to
		self.trans = transformation
		self.refNum = None
class Node:
	'''a node represents a set of observations that
	we think correspond to the same metabolite. E.g. a 
	set of isoforms and adducts'''
	def __init__(self):
		self.name = 'Unnamed Node'
		self.observations = set([])
		self.knowns = set([])
		self.transformations = DSet()
		self.isoforms = DSet()
		self.adducts = DSet()
		# to keep track of when a proposed transformation between two 
		# metabolites is also recorded as an adduct/isoform
		self.nodes_to = []
		self.nodes_from = []
		self.edges_to_self = []
		self.edges_from_self = []
		self.refNum = None
		self.refTransformations = set([])
		self.refAdducts = set([])
	def addAdduct(self, adduct):
		self.adducts.add(adduct)
		self.refAdducts.add(adduct.adductTrans)
	def addIsoform(self, isoform):
		self.isoforms.add(isoform)
	def addTransformation(self, trans):
		self.transformations.add(trans)
		self.refTransformations.add(trans.trans)
	def addKnown(self, known):
		self.knowns.add(known)
	def addObservation(self, obs):
		self.observations.add(obs)
	def merge(node1, node2):
		node_out = Node()
		for var, value in vars(node1).items():
			if type(value) == list:
				setattr(node_out, var, getattr(node1, var) + getattr(node2, var))
			elif var == 'name' or var == 'refNum':
				continue
			else:
				setattr(node_out, var, getattr(node1, var).union(getattr(node2, var)))
		return node_out
	def fromObs(obs):
		n_out = Node()
		n_out.observations.add(obs)
		return n_out
	def addEdgeFrom(self, edge):
		self.nodes_to.append(edge.node_to)
		self.edges_from_self.append(edge)
	def addEdgeTo(self, edge):
		self.nodes_from.append(edge.node_from)
		self.edges_to_self.append(edge)

	# DB must be already connected
	# NOTE: does not include edges/nodes to/from for now
	def fromSqlRow(refNum, db_connected):
		db = db_connected
		n = Node()
		row = [x for x in db.query('select * from nodes where refNum = ' + str(refNum))][0]
		n.refNum = row['refNum']
		n.id = row['id']
		for o_refNum in row['observations'].split(';'):
			obs = Obs.fromSqlRow(o_refNum, db)
			n.addObservation(obs)
		for k_refNum in row['knowns'].split(';'):
			known = Obs.fromSqlRow(k_refNum, db)
			n.addKnown(known)
		for t_refNum in row['transformations'].split(';'):
			trans = Trans.fromSqlRow(t_refNum, db)
			n.addTransformation(trans)
			print(trans)
		for i_refNum in row['isoforms'].split(';'):
			isoform = Isoform.fromSqlRow(i_refNum, db)
			n.addIsoform(isoform)
		for a_refNum in row['adducts'].split(';'):
			adduct = Adduct.fromSqlRow(a_refNum, db)
			n.addAdduct(adduct)


		# obs_q = [Obs.fromSqlRow(x) for x in db.query('select * from observations where ')]



class MergedSet:
	"To keep track of all the unique nodes in our data set"

	# Initializes nodes according to the isoforms and adducts
	def __init__(self):
		# self.edges = []
		self.nodes = []

	def __str__(self):
		for n_i, n in enumerate(self.nodes):
			n.name = 'Node ' + str(n_i)
		s_out = ''
		for n_i, n in enumerate(self.nodes):
			s_out += '------------\n'
			s_out += n.name + '\n'
			s_out += '------------\n'
			s = '\tMetabolites - '
			for o in n.observations:
				s += str(o) + '; '
			if s[-2] == ';':
				s = s[ : -2] + '\n'
			else: 
				s += '\n'
			s += '\tNodes To - '
			for e in n.edges_from_self:
				s += e.node_to.name + '(' + e.trans.trans.name + ')' + '; '
			if s[-2] == ';':
				s = s[ : -2] + '\n'
			else:
				s = s + '\n'
			s += '\tNodes From - '
			for e in n.edges_to_self:
				s += e.node_from.name + '(' + e.trans.trans.name + ')' + '; '
			if s[-2] == ';':	
				s = s[ : -2] + '\n\n'
			else:
				s = s + '\n\n'
			s_out += s
		return s_out

	def addEdge(self, edge):
		edge.node_from.addEdgeFrom(edge)
		edge.node_to.addEdgeTo(edge)
		# self.edges.append(edge)
	def addNode(self, node):
		self.nodes.append(node)
		# self.observations = self.observations.union(node.observations)
		# self.isoforms = self.isoforms.union(node.isoforms.data)
		# self.transformations = self.transformations.union(node.transformations.data)
		# self.adducts = self.adducts.union(node.adducts.data)
		# self.knowns = self.knowns.union(node.knowns)
		# self.edges = self.edges + node.edges_to_self + node.edges_from_self
	def removeNodes(self, r_nodes):
		kept_nodes = []
		for node in self.nodes:
			keep = True
			for r_node in r_nodes:
				if r_node == node:
					keep = False
					break
			if keep:
				kept_nodes.append(node)
		self.nodes = kept_nodes
	def mergeNodes(self, node1, node2):
		newNode = Node.merge(node1, node2)
		self.removeNodes([node1, node2])
		self.addNode(newNode)
		return newNode
	def fromObservations(observations):
		ms_out = MergedSet()
		for obs in observations:
			n = Node.fromObs(obs)
			ms_out.addNode(n)
		return ms_out

	# Gets the from and to nodes
	# corresponding to a given adduct/isoform/transformation
	def getNodes(self, d):
		n_from = None
		n_to = None
		for n in self.nodes:
			n_from_hits = [x for x in n.observations if x.equals(d.obs_from)]
			n_to_hits = [x for x in n.observations if x.equals(d.obs_to)]
			if len(n_from_hits) > 0:
				resolved = True
				if n_from != None or len(n_from_hits) > 1:
					resolved = False
					if n_from != None and np.any([x.equals(y) for (x,y) in zip(n_from.observations, n.observations)]):
							n = self.mergeNodes(n, n_from)
							resolved = True
					elif len(n_from_hits) > 1:
						for i in range(len(n_from_hits)):
							if len(n_from_hits) == 1:
								break
							o_1 = n_from_hits[i]
							for j in range(i, len(n_from_hits)):
								if len(n_from_hits) == 1:
									break
								o_2 = n_from_hits[j]
								if o_1.equals(o_2):
									n.observations.remove(o_2)
									n_from_hits = [x for x in n.observations if x.equals(d.obs_from)]
									n_to_hits = [x for x in n.observations if x.equals(d.obs_to)]
									resolved = True
				if not resolved:
					sys.exit('Error: Multiple node hits found for ' + str(d) + ' obs_from')
				n_from = n
			if len(n_to_hits) > 0:
				resolved = True
				if n_to != None or len(n_to_hits) > 1:
					resolved = False
					if n_to != None and np.any([x.equals(y) for (x,y) in zip(n_to.observations, n.observations)]):
						n = self.mergeNodes(n, n_to)
						resolved = True
					elif len(n_to_hits) > 1:
						for i in range(len(n_to_hits)):
							if len(n_to_hits) == 1:
								break
							o_1 = n_to_hits[i]
							for j in range(i, len(n_to_hits)):
								if len(n_to_hits) == 1:
									break
								o_2 = n_to_hits[j]
								if o_1.equals(o_2):
									print('')
									print('')
									print(o_2 in n.observations)
									print('')
									print('')
									n.observations.remove(o_2)
									n_from_hits = [x for x in n.observations if x.equals(d.obs_from)]
									n_to_hits = [x for x in n.observations if x.equals(d.obs_to)]
									resolved = True
				if not resolved:
					sys.exit('Error: Multiple node hits found for ' + str(d) + ' obs_to')
				n_to = n
		if n_from == None or n_to == None:
			sys.exit('Error: No node hits found for ' + str(d) + ' ' + d.name)
		return n_from, n_to
	def addAdduct(self, adduct):
		n_from, n_to = self.getNodes(adduct)
		if n_from == None or n_to == None:
			sys.exit('Error: Trying to add in adduct with one or more unmatched observations')
		if n_from == n_to:
			n = n_from
		else:
			n = self.mergeNodes(n_from, n_to)
		n.addAdduct(adduct)
	def addAdducts(self, adducts):
		for a in adducts:
			self.addAdduct(a)
	def addIsoform(self, isoform):
		n_from, n_to = self.getNodes(isoform)
		if n_from == None or n_to == None:
			sys.exit('Error: Trying to add in isoform with one or more unmatched observations')
		if n_from == n_to:
			n = n_from
		else:
			n = self.mergeNodes(n_from, n_to)
		n.addIsoform(isoform)
	def addIsoforms(self, isoforms):
		for i in isoforms:
			self.addIsoform(i)
	def addTransformation(self, trans):
		n_from, n_to = self.getNodes(trans)
		if n_from == None or n_to == None:
			sys.exit('Error: Trying to add in transformation with one or more unmatched observations')
		if n_from == n_to:
			n = n_from
			n.addTransformation(trans)
		else:
			newEdge = Edge(n_from, n_to, trans)
			self.addEdge(newEdge)
			n_from.addTransformation(trans)
			n_to.addTransformation(trans)
		if n_from == None or n_to == None:
			sys.exit('Error (lib.py): One of n_from or n_to are none')
	def addTransformations(self, transformations):
		for t in transformations:
			self.addTransformation(t)
	def addKnowns(self, knowns):
		for k in knowns:
			o_match = None
			n_match = None
			for n in self.nodes:
				for o in n.observations:
					if o.equals(k):
						if o_match != None or n_match != None:
							sys.exit('Error: multiple node hits when trying to insert known ' + str(k))
						o_match = o
						n_match = n
						n.addKnown(k)
			if o_match == None or n_match == None:
				sys.exit('Error: no node hits when trying to insert known ' + str(k))
	def getObservations(self):
		observations = []
		for n in self.nodes:
			for o in n.observations:
				observations.append(o)
		return observations
	def getIsoforms(self):
		isoforms = []
		for n in self.nodes:
			for i in n.isoforms.data:
				isoforms.append(i)
		return isoforms
	def getAdducts(self):
		adducts = []
		for n in self.nodes:
			for a in n.adducts.data:
				adducts.append(a)
		return adducts
	def getTransformations(self):
		transformations = []
		for n in self.nodes:
			for t in n.transformations.data:
				# To not overcount transformations
				if t.obs_from.isIn(n.observations):
					transformations.append(t)
		return transformations
	def getKnowns(self):
		knowns = []
		for n in self.nodes:
			for k in n.knowns:
				knowns.append(k)
		return knowns
	def getEdges(self):
		edges = []
		for n in self.nodes:
			for e in n.edges_from_self:
				edges.append(e)
		return edges
	def getRefTransformations(self):
		refTransformations = []
		for n in self.nodes:
			for t in n.refTransformations:
				if t not in refTransformations:
					refTransformations.append(t)
		return refTransformations
	def getRefAdducts(self):
		refAdducts = []
		for n in self.nodes:
			for a in n.refAdducts:
				if a not in refAdducts:
					refAdducts.append(a)
		return refAdducts
	def copy(self):
		c = MergedSet()
		var_dict = vars(self)
		for v in var_dict.keys():
			setattr(c, v, copy.deepcopy(var_dict[v]))
		return c
	def save(self, out_name):
		pkl.dump(self, open(out_name, 'wb'))

	# Writes log file
	def writeLog(self, out_name, args):
		with open(out_name, 'w') as f:
			print('Input DB / table - %s, %s' % (args.knowns_db, args.input_table_name), file = f)
			print('Nodes - ' + str(len(self.nodes)), file = f)
			print('Edges - ' + str(len(self.getEdges())), file = f)
			print('Observations - ' + str(len(self.getObservations())), file = f)
			print('Transformations -' + str(len(self.getTransformations())), file = f)
			print('Knowns - ' + str(len(self.getKnowns())), file = f)
			print('Adducts - ' + str(len(self.getAdducts())), file = f)
			print('Isoforms - ' + str(len(self.getIsoforms())), file = f)
			if args.t_end_time != None:
				print('Merge Runtime - ' + str(round(float(time.time() - args.t_end_time), 2)) + 's', file = f)

	def toSql(self, fname, overwrite = True):
		fname = os.path.abspath(fname)
		if overwrite:
			if os.path.exists(fname):
				os.system('rm ' + fname)
		db = dataset.connect('sqlite:///' + fname)

		observations = self.getObservations()
		adducts = self.getAdducts()
		isoforms = self.getIsoforms()
		transformations = self.getTransformations()
		nodes = self.nodes
		edges = self.getEdges()
		knowns = self.getKnowns()
		refTransformations = self.getRefTransformations()
		refAdducts = self.getRefAdducts()

		all_data = {
		'observations': observations,
			'adducts': adducts,
			'isoforms': isoforms,
			'transformations': transformations,
			'nodes': nodes,
			'edges': edges,
			'knowns' : knowns,
			'refTransformations' : refTransformations,
			'refAdducts' : refAdducts
			}

		prefixes = {
		'observations': 'o',
		'adducts': 'a',
		'isoforms': 'i',
		'transformations': 't',
		'nodes': 'n',
		'edges': 'e',
		'knowns' : 'k',
		'refTransformations' : 'rT',
		'refAdducts' : 'rA'
		}

		column_names = {x : set([]) for x in all_data.keys()}

		# Label with ID numbers
		print('Labeling reference nubmers for DB...')
		for table_name, data_list in all_data.items():
			id_num = 0
			for d in data_list:
				d.refNum = prefixes[table_name] + str(id_num)
				id_num += 1
				column_names[table_name] = column_names[table_name].union(set(vars(d).keys()))
		# Convert to Sql rows and insert into table
		print('Converting objects to sql rows...')
		for table_name, data_list in all_data.items():
			sql_rows = []
			for x in data_list:
				r = toSqlRow(x)
				for c in column_names[table_name]:
					if c not in r.keys():
						r[c] = None
				if 'id' in r:
					del r['id']
				sql_rows.append(r)
			print('Inserting %s' % table_name)
			db[table_name].insert_many(sql_rows, ensure='pid')


		# edge_rows = [x for x in db.query('select * from edges')]
		# edges = []
		# for edge_row in edge_rows:
		# 	[n_from] = [x for x in nodes if str(x.refNum) == edge_row['node_from']]
		# 	[n_to] = [x for x in nodes if str(x.refNum) == edge_row['node_to']]
		# 	if str(edge_row['trans']) not in [str(x.refNum) for x in transformations]:
		# 		continue
		# 	[trans] = [x for x in transformations if str(x.refNum) == edge_row['trans']]

		# 	e_refNum = edge_row['refNum']
		# 	e = Edge(n_from, n_to, trans)
		# 	e.refNum = e_refNum
		# 	edges.append(e)

		# for n in nodes:
		# 	n.nodes_to = [x for x in nodes if str(x.refNum) in n.nodes_to]
		# 	n.nodes_from = [x for x in nodes if str(x.refNum) in n.nodes_from]
		# 	n.edges_to_self = [x for x in edges if str(x.refNum) in n.edges_to_self]
		# 	n.edges_from_self = [x for x in edges if str(x.refNum) in n.edges_from_self]



		# ms_out = MergedSet()
		# ms_out.nodes = nodes
		# return ms_out





	def onlyKnowns(self):
		ms = self.copy()
		for node in ms.nodes:
			keep_node = True
			if not keep_node:
				ms.removeNode(node)
		return ms



################################################################################
# Molecular formula to vector
################################################################################
def atomVector(smiles):
    m = Chem.MolFromSmiles(smiles)
    m = Chem.AddHs(m)
    atoms = m.GetAtoms()
    atom_vec = {}
    for a in atoms:
        if a.GetSymbol() in atom_vec:
            atom_vec[a.GetSymbol()] += 1
        else:
            atom_vec[a.GetSymbol()] = 1
    return atom_vec

def atomsFromObs(obs):
	smiles = obs['SMILES']
	vec = atomVector(smiles)
	return vec

def atomsFromRefTrans(rt):
	atoms = yaml.load(rt['atoms'])
	return atoms


# Gives the vector v2 - v1. NOTE not the same as vecDifference
def vecSubtract(v2, v1):
	out_vec = {x:v2[x] for x in v2}
	for a in v1.keys():
		if a not in v2.keys():
			if v1[a] == 0:
				out_vec[a] = 0
			elif v1[a] < 0:
				out_vec[a] = abs(v1[a])
			else:
				print('v2 - %s' %str(v2))
				print('v1 - %s' %str(v1))
				sys.exit('Error (lib.py): Trying to subtract atoms that dont exist')
		out_vec[a] -= v1[a]
	return out_vec

def vecDifference(v_from, v_to):
	atoms_from = v_from.keys()
	atoms_to = v_to.keys()

	out_vec = {}
	for a in set(atoms_from).symmetric_difference(set(atoms_to)):
		out_vec[a] = 0

	for a in set(atoms_from).intersection(set(atoms_to)):
		out_vec[a] = v_to[a] - v_from[a]

	return out_vec

def vecAdd(v1, v2):
	out_vec = {x:v1[x] for x in v1}
	for a in v2.keys():
		if a in out_vec.keys():
			out_vec[a] += v2[a]
		else:
			out_vec[a] = v2[a]
	return out_vec


def vecGuess(unknown, known, trans, rt):
	k_num = known['refNum']
	trans_atoms = atomsFromRefTrans(rt)
	# known is the obs_from
	if trans['obs_from'] == k_num:
		from_atoms = atomsFromObs(known)
		guess = vecAdd(from_atoms, trans_atoms)
	# known is the obs_to
	else:
		to_atoms = atomsFromObs(known)
		guess = vecSubtract(to_atoms, trans_atoms)

	for a in guess:
		if guess[a] < 0:
			print('\n')
			print('unknown (%s) - %s' % (unknown['refNum'], str(unknown['name'])))
			print('\n\nknown (%s) - %s' % (known['refNum'], str(atomsFromObs(known))))
			print('\n\ntrans - %s' % trans['name'])
			print('\n\nrefTrans - %s' % str(atomsFromRefTrans(rt)))
			print('\n\nguess - %s' % str(guess))
			sys.exit('\nError (lib.py) - Impossible guess')

	return guess

		
def vecEquals(v1, v2):
	atoms_1 = v1.keys()
	atoms_2 = v2.keys()

	for a in set(atoms_1).difference(set(atoms_2)):
		if v1[a] > 0:
			return False
	for a in set(atoms_2).difference(set(atoms_1)):
		if v2[a] > 0:
			return False
	for a in set(atoms_1).intersection(set(atoms_2)):
		if v1[a] != v2[a]:
			return False
	return True









