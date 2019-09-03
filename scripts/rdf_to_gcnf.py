import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from os.path import join, relpath
from glob import glob
import rdflib
import numpy as np

from rdflib import Namespace
from scripts.omics_to_rdf import *
from joblib import delayed
from line_profiler import LineProfiler
from distutils.util import strtobool

class RdfToGcnf(OmicsToRDF):
	def __init__(self, cell_id):
		self.__gdsc_g = None
		self.__ccle_g = None
		self.__cell_id = cell_id
		self.__gdsc_file_path = get_param(DbName.GCN, 'in_gdsc_folder')[0] + GDSC_PREFIX + cell_id + EXT
		self.__ccle_file_path = get_param(DbName.GCN, 'in_ccle_folder')[0] + CCLE_PREFIX + cell_id + EXT
		self.__result_path = get_param(DbName.GCN, 'out_folder')[0]

		self.__variants_ns = Namespace("https://portals.broadinstitute.org/links/variants/")
		self.__faldo_ns = Namespace("http://biohackathon.org/resource/faldo#")
		self.__skos_ns = Namespace("http://www.w3.org/2004/02/skos/core#")

	@property
	def gdsc_graph(self):
		if self.__gdsc_g is None:
			self.__gdsc_g = rdflib.Graph().parse(self.__gdsc_file_path, format='n3')
		return self.__gdsc_g

	@gdsc_graph.setter
	def gdsc_graph(self, value):
		self.__gdsc_g = value

	@property
	def ccle_graph(self):
		if self.__ccle_g is None:
			self.__ccle_g = rdflib.Graph().parse(self.__ccle_file_path, format='n3')
		return self.__ccle_g

	@ccle_graph.setter
	def ccle_graph(self, value):
		self.__ccle_g = value

	def get_query(self, db_name, query_context):
		qres = ''
		if db_name == DbName.CCLE:
			qres = self.ccle_graph.query(query_context)
		elif db_name == DbName.GDSC:
			qres = self.gdsc_graph.query(query_context)
		return qres


	def get_cell_line_info(self):
		# 細胞株情報を取得
		query = """
			PREFIX ontology: <https://www.cancerrxgene.org/translation/gdsc/ontology#>
			PREFIX cclec: <https://portals.broadinstitute.org/links/cell_line/>
			PREFIX m2r: <http://med2rdf/org/ontology/med2rdf#>
			SELECT DISTINCT ?cell_line_name ?cls ?tissue ?tissue_sub
			WHERE {
			  ?cell_line_name rdf:type m2r:CellLine; ontology:tcga_classification ?cls; ontology:tissue ?tissue; ontology:tissue_subtype ?tissue_sub .
			}
			"""
		qres = self.get_query(DbName.GDSC, query)
		if qres == '':
			return 'INVALID'

		cell_stmt = ''
		cls_tis_stmt = ''

		for row in qres:
			cell_stmt = row['cell_line_name'].split('/')[-1]
			cls_stmt = 'TCGAclassification:{}'.format(row['cls'].split('/')[-1].upper())
			tissue_stmt = 'Tissue:{}'.format(row['tissue'].split('/')[-1].lower())
			tissue_sub_stmt = 'TissueSubtype:{}'.format(row['tissue_sub'].split('/')[-1].lower())
			cls_tis_stmt = '\t'.join([cls_stmt, tissue_stmt, tissue_sub_stmt])

		return cell_stmt, cls_tis_stmt

	def get_drug_list(self):
		# 薬剤リストを取得
		query = '''
			PREFIX m2r: <http://med2rdf/org/ontology/med2rdf#>
			SELECT DISTINCT ?drug_id
			WHERE {{
			  ?s1 rdf:type m2r:CellLine; m2r:drug ?o1 .
			  ?o1 m2r:drug ?drug_id
			}}'''
		gdsc_qres = self.get_query(DbName.GDSC, query)
		if gdsc_qres == '':
			return 'INVALID'

		return [drug_num.split('/')[-1] for drug_num in [node['drug_id'] for node in gdsc_qres]]

	def get_drug_info(self, drug_num):
		# 薬剤情報を取得
		query = '''
			PREFIX m2r: <http://med2rdf/org/ontology/med2rdf#>
			PREFIX gdsc: <https://www.cancerrxgene.org/translation/gdsc/>
			PREFIX gdscd: <https://www.cancerrxgene.org/translation/links/drug/>
			PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
			PREFIX sio: <http://semanticscience.org/resource/>
			PREFIX uo: <http://purl.obolibrary.org/obo/>
			
			SELECT DISTINCT ?drug_name ?id ?pathway ?ic50 ?mconc ?auc
			WHERE {{
			  gdscd:{} rdfs:label ?drug_name; skos:altLabel ?id; gdsc:target_pathway ?pathway; sio:SIO_000216 ?o1; sio:SIO_000216 ?o2; sio:SIO_000216 ?o3 .
			  ?o1 rdf:type uo:MI_0641; sio:SIO_000300 ?ic50 .
			  ?o2 rdf:type uo:ScreeningMaxConcentration; sio:SIO_000300 ?mconc .
			  ?o3 rdf:type uo:MCIT_C64774; sio:SIO_000300 ?auc .
			}}'''.format(drug_num)

		qres = self.get_query(DbName.GDSC, query)
		if qres == '':
			return 'INVALID'

		ic50_low_th = float(get_param(DbName.GCN, 'ic50_low_th')[0])
		auc_high_th = float(get_param(DbName.GCN, 'auc_high_th')[0])
		auc_low_th = float(get_param(DbName.GCN, 'auc_low_th')[0])
		id_stmt, name_stmt, pathway_stmt, ic50_stmt, auc_stmt, stmt = '', '', '', '', '', ''
		for row in qres:
			name_stmt = row["drug_name"].lower()
			id_stmt = row["id"]
			pathway_stmt = 'TargetPathway:{}'.format(row["pathway"])

			if ic50_low_th < float(row["ic50"]):
				label_HL = 'IC50:LOW'
			else:
				label_HL = 'IC50:HIGH'

			micro_ic = np.exp(float(row["ic50"]))
			if micro_ic < float(row['mconc']):
				label_SR = 'IC50:sensitive'
			else:
				label_SR = 'IC50:resistant'

			ic50_stmt = '\t'.join([label_HL, label_SR])

			if float(row["auc"]) < auc_high_th:
				auc_stmt = 'AUC:HIGH'
			elif (auc_high_th <= float(row["auc"])) & (float(row["auc"]) <= auc_low_th):
				auc_stmt = 'AUC:MIDDLE'
			else:
				auc_stmt = 'AUC:LOW'
		stmt = '\t'.join([id_stmt, name_stmt, pathway_stmt, ic50_stmt, auc_stmt])
		return stmt

	def get_gene_list(self):
		# 遺伝子リストを取得
		query = '''
			PREFIX m2r: <http://med2rdf/org/ontology/med2rdf#>
			SELECT DISTINCT ?ensembl ?symbol
			WHERE {{
			  ?s1 rdf:type m2r:CellLine ; m2r:gene ?o .
			  ?o rdfs:seeAlso ?ensembl; rdfs:label ?symbol
			}}'''

		ccle_qres = self.get_query(DbName.CCLE, query)
		if ccle_qres == '':
			return 'INVALID'

		ensembl_list = [gene_num.split('/')[-1] for gene_num in [node['ensembl'] for node in ccle_qres]]
		symbol_list = [gene_num.split('/')[-1] for gene_num in [node['symbol'] for node in ccle_qres]]

		return ensembl_list, symbol_list

	def get_variants_list(self, ensmbl):
		# 変異リストを取得
		query = '''
			PREFIX m2r: <http://med2rdf/org/ontology/med2rdf#>
			PREFIX faldo: <http://biohackathon.org/resource/faldo#>
			SELECT DISTINCT ?chr ?ref ?alt ?start ?vtype ?vclass
			WHERE {{
			  ?s1 rdfs:seeAlso "{}"; m2r:variation ?o1 .
			  ?o1 faldo:location ?o2; m2r:reference_allele ?ref; m2r:alternative_allele ?alt; ontology:variant_type ?vtype; ontology:variant_classification ?vclass .
			  ?o2 faldo:begin ?o3 .
			  ?o3 faldo:reference ?chr; faldo:position ?start .
			}}'''.format(ensmbl)

		qres = self.get_query(DbName.CCLE, query)

		stmt_list = []
		for row in qres:
			chr = row["chr"].split('/')[-1]
			ref = row["ref"].split('/')[-1]
			alt = row["alt"].split('/')[-1]
			start = row["start"].split('/')[-1]
			vtype = row["vtype"].split('/')[-1]
			vcls = row["vclass"].split('/')[-1]
			stmt_list.append('\t'.join(["{}_{}_{}_{}".format(chr, start, ref, alt),
			                          'VariantClassification:{}'.format(vcls.lower()),
			                          'VariantType:{}'.format(vtype.upper())]))
		return stmt_list

	def create_gcn_file(self):
		output_path = self.__result_path + 'result_gcn_{}.txt'.format(self.__cell_id)
		if os.path.isfile(output_path):
			return

		cell_line_stmt, tissue_stmt = self.get_cell_line_info()
		drug_list = self.get_drug_list()
		ensembl_list, symbol_list = self.get_gene_list()

		for drug in drug_list:
			drug_stmt = self.get_drug_info(drug)
			all_stmt = []
			for ensembl, symbol in zip(ensembl_list, symbol_list):
				ensembl_stmt = "ENSEMBL:" + ensembl
				symbol_stmt = "{}_HUMAN".format(symbol)
				variants_list = self.get_variants_list(ensembl)
				for variants_stmt in variants_list:
					stmt = [cell_line_stmt, variants_stmt, symbol_stmt, ensembl_stmt, tissue_stmt, drug_stmt]
					all_stmt.append('\t'.join(stmt))

			if len(all_stmt) > 0:
				result = '\n'.join(all_stmt)
				result += '\n'
				with open(output_path, 'a') as f:
					f.write(result)

				del all_stmt
				del result


def process(cell_id):
	if strtobool(get_param(DbName.COMMON, 'profile')[0]):
		print('cell id : {}'.format(cell_id))
		pr = LineProfiler()
		rtg = RdfToGcnf(cell_id)
		pr.add_function(rtg.create_gcn_file)
		pr.enable()
		rtg.create_gcn_file()
		pr.disable()
		pr.print_stats()
		del rtg
	else:
		rtg = RdfToGcnf(cell_id)
		rtg.create_gcn_file()
		del rtg


if __name__ == '__main__':
	ccle_path = get_param(DbName.GCN, 'in_ccle_folder')[0]
	gdsc_path = get_param(DbName.GCN, 'in_gdsc_folder')[0]

	# ccleの細胞株とgdscの細胞株で共通の細胞株を取得する
	ccle_list = [relpath(x, ccle_path)[len(CCLE_PREFIX):-len(EXT)] for x in glob(join(ccle_path, '*'))]
	gdsc_list = [relpath(x, gdsc_path)[len(GDSC_PREFIX):-len(EXT)] for x in glob(join(gdsc_path, '*'))]
	common_list = list(set(ccle_list) & set(gdsc_list))
	jobs = int(get_param(DbName.COMMON, 'n_jobs')[0])
	print("CCLEとGDSCのRDF/turtleファイルからGCNファイルに変換")
	aprun = ParallelExecutor(n_jobs=jobs)
	aprun(total=len(common_list))(delayed(process)(cell_id) for cell_id in common_list)

