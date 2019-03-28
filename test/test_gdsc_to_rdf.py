from unittest import TestCase
from scripts.gdsc_to_rdf import GdscToRDF
from scripts.omics_to_rdf import *
from rdflib import Literal, Namespace
from SPARQLWrapper import SPARQLWrapper
import pickle

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

ic50_tag_list = ["Cell line name",
		                        "Drug name",
		                        "Drug Id",
		                        "IC50",
		                        "Cosmic sample Id",
		                        "TCGA classification",
		                        "Tissue",
		                        "Tissue sub-type",
		                        "AUC",
		                        "Dataset version",
		                        "IC Result ID"]



class TestGdscToRDF(TestCase):

	def setUp(self):
		with open(get_param(DbName.GDSC, 'in_gdsc_drug_file')[0], 'rb') as f:
			self.ic50 = pickle.load(f)

		self.output_path = get_param(DbName.GDSC, 'out_folder')[0] + 'debug_rdf.txt'

		self.__gdsc_ns = Namespace("https://www.cancerrxgene.org/translation/Drug/")
		self.__skos_ns = Namespace("http://www.w3.org/2004/02/skos/core#")
		self.sparql = SPARQLWrapper(endpoint='http://sparql.hegroup.org/sparql/', returnFormat='tsv')

		self.graph = None

	def tearDown(self):
		if self.graph is not None:
			self.graph.serialize(destination=self.output_path, format="turtle")
			with open(self.output_path) as f:
				result = f.read()
			print(result)

	def test_should_create_DMS53_cell_line_graph(self):
		row = self.ic50[self.ic50['Cosmic sample Id'] == 907295]
		gtr = GdscToRDF(row)
		self.graph = gtr.create_gdsc_single_cell_turtle(True)

		expected = [Literal("dms53_doxorubicin"), Literal("dms53_etoposide"), Literal("dms53_gemcitabine")]
		result = []
		for s1, p1, o1 in self.graph.triples((None, self.__skos_ns["altLabel"], None)):
			if ("_doxorubicin" in o1) | ("_etoposide" in o1) | ("_gemcitabine" in o1):
				result.append(o1)
		result.sort()
		expected.sort()

		self.assertListEqual(result, expected)

	def test_should_create_22rv1_cell_line_graph(self):
		row = self.ic50[self.ic50['Cosmic sample Id'] == 924100]
		gtr = GdscToRDF(row)
		self.graph = gtr.create_gdsc_single_cell_turtle(True)

		expected = [Literal("22rv1_doxorubicin"), Literal("22rv1_etoposide"), Literal("22rv1_gemcitabine")]
		result = []
		for s1, p1, o1 in self.graph.triples((None, self.__skos_ns["altLabel"], None)):
			if ("_doxorubicin" in o1) | ("_etoposide" in o1) | ("_gemcitabine" in o1):
				result.append(o1)
		result.sort()
		expected.sort()

		self.assertListEqual(result, expected)

	def test_should_be_equal_to_auc_value_of_DMS53_cell_line(self):
		row = self.ic50[self.ic50['Cosmic sample Id'] == 907295]
		gtr = GdscToRDF(row)
		self.graph = gtr.create_gdsc_single_cell_turtle(True)

		query_context = '''
			PREFIX gdscd: <https://www.cancerrxgene.org/translation/links/drug/>
			PREFIX sio: <http://semanticscience.org/resource/>
			PREFIX uo: <http://purl.obolibrary.org/obo/>
			SELECT DISTINCT ?auc
			WHERE {{
			  gdscd:{} sio:SIO_000216 ?o1 .
			  ?o1 rdf:type uo:MCIT_C64774; sio:SIO_000300 ?auc .
			}}'''.format('169396')

		qres = self.graph.query(query_context)

		for row in qres:
			result = float(row["auc"])

		expected = 0.89359217
		self.assertEqual(result, expected)

	def test_should_be_equal_to_ic50_value_of_DMS53_cell_line(self):
		row = self.ic50[self.ic50['Cosmic sample Id'] == 907295]
		gtr = GdscToRDF(row)
		self.graph = gtr.create_gdsc_single_cell_turtle(True)

		query_context = '''
			PREFIX gdscd: <https://www.cancerrxgene.org/translation/links/drug/>
			PREFIX sio: <http://semanticscience.org/resource/>
			PREFIX uo: <http://purl.obolibrary.org/obo/>
			SELECT DISTINCT ?ic50
			WHERE {{
			  gdscd:{} sio:SIO_000216 ?o1 .
			  ?o1 rdf:type uo:MI_0641; sio:SIO_000300 ?ic50 .
			}}'''.format('169396')

		qres = self.graph.query(query_context)

		for row in qres:
			result = float(row["ic50"])

		expected = 0.39437767
		self.assertEqual(result, expected)

	def test_should_be_equal_to_auc_value_of_22rv1_cell_line(self):
		row = self.ic50[self.ic50['Cosmic sample Id'] == 924100]
		gtr = GdscToRDF(row)
		self.graph = gtr.create_gdsc_single_cell_turtle(True)

		query_context = '''
			PREFIX gdscd: <https://www.cancerrxgene.org/translation/links/drug/>
			PREFIX sio: <http://semanticscience.org/resource/>
			PREFIX uo: <http://purl.obolibrary.org/obo/>
			SELECT DISTINCT ?auc
			WHERE {{
			  gdscd:{} sio:SIO_000216 ?o1 .
			  ?o1 rdf:type uo:MCIT_C64774; sio:SIO_000300 ?auc .
			}}'''.format('299906')

		qres = self.graph.query(query_context)

		for row in qres:
			result = float(row["auc"])

		expected = 0.53709514
		self.assertEqual(result, expected)

	def test_should_be_equal_to_ic50_value_of_22rv1_cell_line(self):
		row = self.ic50[self.ic50['Cosmic sample Id'] == 924100]
		gtr = GdscToRDF(row)
		self.graph = gtr.create_gdsc_single_cell_turtle(True)

		query_context = '''
			PREFIX gdscd: <https://www.cancerrxgene.org/translation/links/drug/>
			PREFIX sio: <http://semanticscience.org/resource/>
			PREFIX uo: <http://purl.obolibrary.org/obo/>
			SELECT DISTINCT ?ic50
			WHERE {{
			  gdscd:{} sio:SIO_000216 ?o1 .
			  ?o1 rdf:type uo:MI_0641; sio:SIO_000300 ?ic50 .
			}}'''.format('299906')

		qres = self.graph.query(query_context)

		for row in qres:
			result = float(row["ic50"])

		expected = 0.25828248
		self.assertEqual(result, expected)
		
	
