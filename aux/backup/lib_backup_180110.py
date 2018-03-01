################################################################################
# Gabe Reder - gkreder@gmail.com
# 08/2017
################################################################################
import pickle as pkl
import os
import math
import time
import sys
from xlrd import open_workbook
import copy
import dataset
################################################################################
# Objects/Classes
################################################################################
class Obs:
	"Class to define observations from input files"
	# assuming input tsv line
	def __init__(self, line, input_type = 'meyer'):
		if input_type == 'meyer':
			self.mz = float(line.split('\t')[2])
			self.rt = float(line.split('\t')[3])
			self.method = line.split('\t')[1]
			self.metabolite = line.split('\t')[0]
			self.metab = self.metabolite
			self.name = self.metab
			self.refNum = None
			self.known = False
		elif input_type == 'hirschhorn':
			self.mz = float(line.split('\t')[2])
			self.rt = float(line.split('\t')[3])
			self.method = line.split('\t')[1]
			self.metabolite = line.split('\t')[0]
			self.metab = self.metabolite
			self.refNum = None
			self.known = False
		else:
			sys.exit('Error: Unrecognized input observation type')
	def __str__(self):
		return self.name

	# takes in file and returns list of observations in file
	def fromFile(fname, input_type = 'meyer', header = True):
		with open(fname, 'r') as f:
			lines = [x.strip() for x in f]
		if header:
			lines = lines[1 : ]

		obs_out = [Obs(x, input_type=input_type) for x in lines]
		return obs_out

	def toSqlRow(self):
		var_dict = vars(self)
		row = {v : str(var_dict[v]) for v in var_dict if v != 'metabolite' and v != 'metab'}
		return row

	def equals(self, o):
		e = True
		if self.mz != o.mz:
			e = False
		if self.rt != o.rt:
			e = False
		if self.method != o.method:
			e = False
		if self.metabolite != o.metabolite:
			e = False
		if self.metab != o.metab:
			e = False
		if self.name != o.name:
			e = False
		return e

	def isKnown(self):
		return self.known


class RefTrans:
	"Class to define input (reference) transformations"
	def __init__(self, name, mz):
		self.name = name
		self.dmz = float(mz)

	# get transformations from excel file
	def fromExcel(e_file):
		has_isoforms = False
		wb = open_workbook(e_file)
		transformation_sheet = wb.sheet_by_name('Common chemical relationships')
		mz_list = [x.value for x in transformation_sheet.col(7)[2 : ]]
		names = [x.value for x in transformation_sheet.col(0)[2 : ]]
		for n in names:
			if n.lower() == 'isoform':
				has_isoforms = True
		transformations = list(zip(mz_list, names))
		transformations = [RefTrans(n, m) for (m,n) in transformations]
		if not has_isoforms:
			sys.exit('Error: please add isoform to chemical trasnformations sheet')
		return transformations
class Trans:
	"Class to define output (proposed) transformations"
	def __init__(self, closest_match, obs_from, obs_to):
		delta_mz_obs = obs_to.mz - obs_from.mz
		# proposed transformation (Trans object)
		self.trans = closest_match
		# delta_mz observed
		self.dmz_obs = delta_mz_obs
		# metabolite from
		self.obs_from = obs_from
		# metabolite to
		self.obs_to = obs_to
		self.refNum = None
		# method
		self.method = obs_from.method
		if obs_to.method != self.method:
			sys.exit('Error: Different observation methods for proposed transformation')

	# Return tsv line given a Trans object
	def line(self, delimiter = '\t'):
		l = []
		l.append(self.obs_from.metab + '--->' + self.obs_to.metab)
		l.append(self.trans.name)
		l.append(str(self.dmz_obs))
		l.append(str(self.trans.dmz - self.dmz_obs))
		l.append(self.obs_from.metab)
		l.append(str(self.obs_from.mz))
		l.append(str(self.obs_from.rt))
		l.append(self.obs_to.metab)
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
			print('Infile - ' + args.in_file, file = f)
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

	def toSqlRow(self):
		r = {}
		r['transformation'] = self.trans.name
		r['dmz'] = str(self.dmz_obs)
		r['obs_from_refNum'] = str(self.obs_from.refNum)
		r['obs_to_refNum'] = str(self.obs_to.refNum)
		r['obs_from'] = self.obs_from.name
		r['obs_to'] = self.obs_to.name
		r['refNum'] = str(self.refNum)
		r['method'] = self.method
		return r

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
class AdductTrans:
	"Proposed found transformation corresponding to two adducts of a metabolite"
	def __init__(self, closest_match, obs_from, obs_to):
		delta_mz_obs = obs_to.mz - obs_from.mz
		self.adductTrans = closest_match
		self.obs_from = obs_from
		self.obs_to = obs_to
		self.dmz_obs = delta_mz_obs
		self.method = obs_from.method
		if obs_to.method != self.method:
			sys.exit('Error: Different observation methods for proposed adduct transformation')
		self.mode = closest_match.mode
	# Return tsv line given a PropTrans object
	def line(self, delimiter = '\t'):
		l = []
		l.append(self.obs_from.metab + '--->' + self.obs_to.metab)
		l.append(self.adductTrans.name)
		l.append(str(self.dmz_obs))
		l.append(str(self.adductTrans.dmz - self.dmz_obs))
		l.append(self.obs_from.metab)
		l.append(str(self.obs_from.mz))
		l.append(str(self.obs_from.rt))
		l.append(self.obs_to.metab)
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
			print('Infile - ' + args.in_file, file = f)
			print('M/Z Tolerance - ' + str(args.adduct_tol_mz), file = f)
			print('RT Tolerance - ' + str(args.adduct_tol_rt), file = f)
			print('Input observations - ' + str(len(lines)), file = f)
			print('Found Adducts - ' + str(len(saved_adducts)), file = f)

	def toSqlRow(self):
		r = {}
		r['adductTrans'] = str(self.adductTrans.name)
		r['obs_from_refNum'] = str(self.obs_from.refNum)
		r['obs_to_refNum'] = str(self.obs_to.refNum)
		r['obs_from'] = self.obs_from.name
		r['obs_to'] = self.obs_to.name
		r['dmz'] = str(self.dmz_obs)
		r['method'] = str(self.method)
		return r

class Isoform():
	"Isoform meaning two observations that look the same"
	def __init__(self, obs_from, obs_to):
		self.obs_from = obs_from
		self.obs_to = obs_to
		self.dmz_obs = obs_to.mz - obs_from.mz
		self.drt = obs_to.rt - obs_from.rt
		self.method = obs_from.method
		if obs_to.method != self.method:
			sys.exit('Error: Different observation methods for proposed isoform')

	# Return tsv line given a Isoform object
	def line(self, delimiter = '\t'):
		l = []
		l.append(self.obs_from.metab + '--->' + self.obs_to.metab)
		l.append(str(self.dmz_obs))
		l.append(str(0.0 - self.dmz_obs))
		l.append(self.obs_from.metab)
		l.append(str(self.obs_from.mz))
		l.append(str(self.obs_from.rt))
		l.append(self.obs_to.metab)
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
			print('Infile - ' + args.in_file, file = f)
			print('M/Z Tolerance - ' + str(args.adduct_tol_mz), file = f)
			print('RT Tolerance - ' + str(args.adduct_tol_rt), file = f)
			print('Input observations - ' + str(len(lines)), file = f)
			print('Found Isoforms - ' + str(len(saved_isoforms)), file = f)

	def toSqlRow(self):
		r = {}
		r['obs_from_refNum'] = str(self.obs_from.refNum)
		r['obs_to_refNum'] = str(self.obs_to.refNum)
		r['obs_from'] = self.obs_from.name
		r['obs_to'] = self.obs_to.name
		r['dmz'] = str(self.dmz_obs)
		r['drt'] = str(self.drt)
		r['method'] = str(self.method)
		return r

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


class Edge:
	"connecting two nodes"
	def __init__(self, node_from, node_to, transformation):
		self.node_from = node_from
		self.node_to = node_to
		self.trans = transformation
		self.refNum = None

	def toSqlRow(self):
		r = {}
		r['node_from'] = str(self.node_from.refNum)
		r['node_to'] = str(self.node_to.refNum)
		r['trans'] = str(self.trans.refNum)
		r['refNum'] = str(self.refNum)
		return r

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
		self.internal_transformations = DSet()
		self.nodes_to = []
		self.nodes_from = []
		self.edges_to_self = []
		self.edges_from_self = []
		self.refNum = None
	# def __str__(self):
	# 	out_string = ''
	# 	out_string += 'Observations:\n\t'
	# 	for o in self.observations:
	# 		out_string += o.name + ','
	# 	return out_string
	def printObservations(self):
		p_set = {str(o) for o in self.observations}
		print(p_set)
	def addAdduct(self, adduct):
		self.adducts.add(adduct)
		self.observations.add(adduct.obs_from)
		self.observations.add(adduct.obs_to)
	def addIsoform(self, isoform):
		self.isoforms.add(isoform)
		self.observations.add(isoform.obs_from)
		self.observations.add(isoform.obs_to)
	def addInternalTransformation(self, trans):
		self.internal_transformations.add(trans)
	def addTransformation(self, trans):
		self.transformations.add(trans)
	def addKnown(self, known):
		self.knowns.add(known)
	def removeKnown(self, known):
		self.knowns.remove(known)
	def containsFrom(self, transformation):
		obs_from = transformation.obs_from
		if obs_from in self.observations:
			return True
		else:
			return False
	def containsTo(self, transformation):
		obs_to = transformation.obs_to
		if obs_to in self.observations:
			return True
		else:
			return False, 'None'
	def containAdduct(self, adduct):
		if adduct.obs_from in self.observations:
			return True
		if adduct.obs_to in self.observations:
			return True
		else:
			return False
	def merge(node1, node2):
		node_out = Node()
		node_out.observations = node1.observations.union(node2.observations)
		node_out.transformations = node1.transformations.union(node2.transformations)
		node_out.isoforms = node1.isoforms.union(node2.isoforms)
		node_out.adducts = node1.adducts.union(node2.adducts)
		node_out.internal_transformations = node1.internal_transformations.union(node2.internal_transformations)
		node_out.nodes_to = node1.nodes_to + node2.nodes_to
		node_out.nodes_from = node1.nodes_from + node2.nodes_from
		node_out.edges_to_self = node1.edges_to_self + node2.edges_to_self
		node_out.edges_from_self = node1.edges_from_self + node2.edges_from_self
		node_out.knowns = node1.knowns.union(node2.knowns)
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
	def removeObs(self, obs):
		self.observations.remove(obs)

	def toSqlRow(self):
		n_row = {}
		n_row['nodes_from'] = ';'.join([str(x.refNum) for x in self.nodes_from])
		n_row['nodes_to'] = ';'.join([str(x.refNum) for x in self.nodes_to])
		n_row['edges_to_self'] = ';'.join([str(x.refNum) for x in self.edges_to_self])
		n_row['edges_to_self'] = ';'.join([str(x.refNum) for x in self.edges_to_self])
		n_row['observations'] = ';'.join([str(x.refNum) for x in self.observations])
		n_row['adducts'] = ';'.join([str(x.refNum) for x in self.adducts.data])
		n_row['isoform'] = ';'.join([str(x.refNum) for x in self.isoforms.data])
		n_row['transformations'] = ';'.join([str(x.refNum) for x in self.transformations.data])
		n_row['internal_transformations'] = ';'.join([str(x.refNum) for x in self.internal_transformations.data])
		return n_row

class MergedSet:
	"To keep track of all the unique nodes in our data set"

	# Initializes nodes according to the isoforms and adducts
	def __init__(self):
		self.edges = []
		self.nodes = []
		self.observations = set([])
		self.transformations = set([])
		self.adducts = set([])
		self.isoforms = set([])
		self.knowns = set([])

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
		self.edges.append(edge)
	def addNode(self, node):
		self.nodes.append(node)
		self.observations = self.observations.union(node.observations)
		self.isoforms = self.isoforms.union(node.isoforms.data)
		self.transformations = self.transformations.union(node.transformations.data)
		self.adducts = self.adducts.union(node.adducts.data)
		self.knowns = self.knowns.union(node.knowns)
		self.edges = self.edges + node.edges_to_self + node.edges_from_self
	def removeNode(self, r_node):
		nodes_out = []
		for node in self.nodes:
			if r_node == node:
				continue
			else:
				nodes_out.append(node)
		self.observations = self.observations.difference(r_node.observations)
		self.transformations = self.transformations.difference(r_node.transformations.data)
		self.isoforms = self.isoforms.difference(r_node.isoforms.data)
		self.adducts = self.adducts.difference(r_node.adducts.data)
		remove_knowns = []
		for k in self.knowns:
			for y in r_node.knowns:
				if k.equals(y):
					remove_knowns.append(k)
		self.knowns = set([x for x in self.knowns if x not in remove_knowns])
		for e_from in r_node.edges_from_self:
			if e_from in self.edges:
				self.edges.remove(e_from)
		for e_to in r_node.edges_to_self:
			if e_to in self.edges:
				self.edges.remove(e_to)
		self.nodes = nodes_out
	def mergeNodes(self, node1, node2):
		newNode = Node.merge(node1, node2)
		self.removeNode(node1)
		self.removeNode(node2)
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
			if d.obs_from in n.observations:
				n_from = n
			if d.obs_to in n.observations:
				n_to = n
		return n_from, n_to
	def addAdduct(self, adduct):
		self.adducts.add(adduct)
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
		self.isoforms.add(isoform)
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
		self.transformations.add(trans)
		n_from, n_to = self.getNodes(trans)
		if n_from == None or n_to == None:
			sys.exit('Error: Trying to add in transformation with one or more unmatched observations')
		if n_from == n_to:
			n = n_from
			n.addInternalTransformation(trans)
		else:
			newEdge = Edge(n_from, n_to, trans)
			self.addEdge(newEdge)
			n_from.addTransformation(trans)
			n_to.addTransformation(trans)
	def addTransformations(self, transformations):
		for t in transformations:
			self.addTransformation(t)
	def addKnowns(self, knowns):
		for k in knowns:
			o_match = None
			n_match = None
			for counter, o in enumerate(self.observations):
				if k.equals(o):
					if o_match != None:
						sys.exit('Error: Found multiple observation matches when adding known')
					o_match = o
			if o_match == None:
				sys.exit('Error: Found no observation matches when adding known')
			self.knowns.add(o_match)
			o_match.known = True
			for n in self.nodes:
				for o in n.observations:
					if o.equals(o_match):
						n.addKnown(k)
						n_match = n
						break
			if n_match == None:
				sys.exit('Error: Couldnt find a node with a matching observation for this known')

	def countObservations(self):
		numObs = 0
		for node in self.nodes:
			numObs += len(node.observations)
		if numObs != len(self.observations):
			sys.exit('Error: Inconsistent observation set size')
		return len(self.observations)
	
	# return a copy of the data structure
	def copy(self):
		c = MergedSet()
		var_dict = vars(self)
		for v in var_dict.keys():
			setattr(c, v, copy.deepcopy(var_dict[v]))
		return c

	# Saves MergedSet object
	def save(self, out_name):
		pkl.dump(self, open(out_name, 'wb'))

	# Writes log file
	def writeLog(self, out_name, args):
		with open(out_name, 'w') as f:
			print('Input File - ' + args.in_file, file = f)
			print('Nodes - ' + str(len(self.nodes)), file = f)
			print('Edges - ' + str(len(self.edges)), file = f)
			print('Metabolites - ' + str(len(self.observations)), file = f)
			if args.t_end_time != None:
				print('Merge Runtime - ' + str(round(float(time.time() - args.t_end_time), 2)) + 's', file = f)

	def writeSql(self, fname, overwrite = True):
		fname = os.path.abspath(fname)
		if overwrite:
			if os.path.exists(fname):
				os.system('rm ' + fname)

		o_rows = []
		n_rows = []
		a_rows = []
		e_rows = []
		i_rows = []
		t_rows = []
		k_rows = []
		it_rows = []
		i_id = 0
		a_id = 0
		t_id = 0
		it_id = 0
		e_id = 0
		o_id = 0

		# Prime observations, Nodes, edges, isoforms, adducts, by giving them
		# IDs to refer to them by
		for n_id, n in enumerate(self.nodes):
			n.refNum = n_id
			for it in n.internal_transformations.data:
				it.refNum = it_id
				it_id += 1
			for a in n.adducts.data:
				a.refNum = a_id
				a_id += 1
			for i in n.isoforms.data:
				i.refNum = i_id
				i_id += 1
			for t in n.transformations.data:
				t.refNum = t_id
				t_id += 1
			for e in n.edges_to_self:
				e.refNum = e_id
				e_id += 1
			for o in n.observations:
				o.refNum = o_id
				o_id += 1

		db = dataset.connect('sqlite:///' + fname)

		for n in self.nodes:
			n_rows.append(n.toSqlRow())
			for a in n.adducts.data:
				a_rows.append(a.toSqlRow())
			for i in n.isoforms.data:
				i_rows.append(i.toSqlRow())
			for t in n.transformations.data:
				t_rows.append(t.toSqlRow())
			for it in n.internal_transformations.data:
				it_rows.append(it.toSqlRow())
			# Don't double count
			for e in n.edges_to_self:
				e_rows.append(e.toSqlRow())
			for o in n.observations:
				o_row = o.toSqlRow()
				o_rows.append(o_row)
				for k in n.knowns:
					if k.equals(o):
						k_rows.append(o_row)
						break

		# Get rid of repeats
		o_rows_temp = []
		n_rows_temp = []
		a_rows_temp = []
		e_rows_temp = []
		i_rows_temp = []
		t_rows_temp = []
		k_rows_temp = []
		it_rows_temp = []

		for r in o_rows:
			add = True
			for r_temp in o_rows_temp:
				if r == r_temp:
					add = False
					break
			if add:
				o_rows_temp.append(r)
		for r in n_rows:
			add = True
			for r_temp in n_rows_temp:
				if r == r_temp:
					add = False
					break
			if add:
				n_rows_temp.append(r)
		for r in a_rows:
			add = True
			for r_temp in a_rows_temp:
				if r == r_temp:
					add = False
					break
			if add:
				a_rows_temp.append(r)
		for r in e_rows:
			add = True
			for r_temp in e_rows_temp:
				if r == r_temp:
					add = False
					break
			if add:
				e_rows_temp.append(r)
		for r in i_rows:
			add = True
			for r_temp in i_rows_temp:
				if r == r_temp:
					add = False
					break
			if add:
				i_rows_temp.append(r)
		for r in t_rows:
			add = True
			for r_temp in t_rows_temp:
				if r == r_temp:
					add = False
					break
			if add:
				t_rows_temp.append(r)
		for r in k_rows:
			add = True
			for r_temp in k_rows_temp:
				if r == r_temp:
					add = False
					break
			if add:
				k_rows_temp.append(r)
		for r in it_rows:
			add = True
			for r_temp in it_rows_temp:
				if r == r_temp:
					add = False
					break
			if add:
				it_rows_temp.append(r)

		o_rows = o_rows_temp
		n_rows = n_rows_temp
		a_rows = a_rows_temp
		e_rows = e_rows_temp
		i_rows = i_rows_temp
		t_rows = t_rows_temp
		k_rows = k_rows_temp
		it_rows = it_rows_temp

		db['nodes'].insert_many(n_rows)
		db['transformations'].insert_many(t_rows)
		db['adducts'].insert_many(a_rows)
		db['isoforms'].insert_many(i_rows)
		db['knowns'].insert_many(k_rows)
		db['observations'].insert_many(o_rows)
		db['edges'].insert_many(e_rows)
		db['internal_transformations'].insert_many(it_rows)

		


	def onlyKnowns(self):
		ms = self.copy()
		for node in ms.nodes:
			keep_node = True
			if not keep_node:
				ms.removeNode(node)
		return ms


################################################################################
# PRINT TIME ELAPSED
################################################################################
def check_time(i, lines, start_time, tabs = 0):
	index_1 = int(math.ceil((len(lines) / 100)))
	index_10 = int(index_1 * 10)
	index_25 = int(index_1 * 25)
	index_50 = int(index_1 * 50)
	index_75 = int(index_1 * 75)
	index_90 = int(index_1 * 90)
	print_string = ''
	for b in range(tabs):
		print_string += '\t'
	if i == index_1:
		print(print_string + '1% --- ' + 
			  str(round(float(time.time() - start_time),2)) + 's')
		sys.stdout.flush()
	if i == index_10:
		print(print_string + '10% --- ' + 
			  str(round(float(time.time() - start_time),2)) + 's')
		sys.stdout.flush()
	if i == index_25:
		print(print_string + '25% --- ' + 
			  str(round(float(time.time() - start_time),2)) + 's')
		sys.stdout.flush()
	if i == index_50:
		print(print_string + '50% --- ' + 
			  str(round(float(time.time() - start_time),2)) + 's')
		sys.stdout.flush()
	if i == index_75:
		print(print_string + '75% --- ' + 
			  str(round(float(time.time() - start_time),2)) + 's')
		sys.stdout.flush()
	if i == index_90:
		print(print_string + '90% --- ' + 
			  str(round(float(time.time() - start_time),2)) + 's')
		sys.stdout.flush()
################################################################################
################################################################################