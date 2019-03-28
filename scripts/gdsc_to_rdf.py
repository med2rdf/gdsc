import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from rdflib import Literal, Namespace, Graph, BNode
from rdflib.namespace import RDF, RDFS, DCTERMS, XSD
import pickle
from scripts.omics_to_rdf import *
from joblib import delayed
import re

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


class GdscToRDF(OmicsToRDF):
	def __init__(self, ic50_row):
		super().__init__()
		self.__result_path = get_param(DbName.GDSC, 'out_folder')[0]
		self.__gdsc_ns = Namespace("https://www.cancerrxgene.org/translation/gdsc/")
		self.__cellline_ns = Namespace("https://www.cancerrxgene.org/translation/links/cell_line/")
		self.__drug_ns = Namespace("https://www.cancerrxgene.org/translation/links/drug/")
		self.__cellosaurus_ns = Namespace("https://web.expasy.org/cgi-bin/cellosaurus/")
		self.__m2r_ns = Namespace("http://med2rdf/org/ontology/med2rdf#")
		self.__ontology_ns = Namespace("https://www.cancerrxgene.org/translation/gdsc/ontology#")
		self.__sio_ns = Namespace("http://semanticscience.org/resource/")
		self.__uo_ns = Namespace("http://purl.obolibrary.org/obo/")
		self.__tgo_ns = Namespace("http://purl.jp/bio/101/opentggates/ontology/")
		self.__tgexpc_ns = Namespace("http://purl.jp/bio/101/opentggates/ExperimentalCondition/")
		self.__bao_ns = Namespace("http://www.bioassayontology.org/bao#")
		self.__skos_ns = Namespace("http://www.w3.org/2004/02/skos/core#")

		self.__id_map = None
		self.__ic50 = ic50_row

		self.g = None
		self.init_graph()

		if os.path.isfile(self.__result_path):
			os.remove(self.__result_path)

	@property
	def id_map(self):
		if self.__id_map is None:
			with open(get_param(DbName.COMMON, 'in_cellosaurus_gdsc_file')[0], 'rb') as f:
				self.__id_map = pickle.load(f)
		return self.__id_map

	@id_map.setter
	def id_map(self, value):
		self.__id_map = value

	@property
	def ic50(self):
		if self.__ic50 is None:
			with open(get_param(DbName.GDSC, 'in_gdsc_drug_file')[0], 'rb') as f:
				self.__ic50 = pickle.load(f)

			self.__ic50 = self.__ic50.loc[:, ic50_tag_list]
		return self.__ic50

	@ic50.setter
	def ic50(self, value):
		self.__ic50 = value

	def save_turtle(self, name):
		name = re.sub(r'\W', "_", name)
		out_path = self.__result_path + GDSC_PREFIX + name + EXT
		idx = 2
		while(os.path.isfile(out_path)):
			out_path = self.__result_path + GDSC_PREFIX + name + str(idx) + EXT
			print(out_path)
			idx += 1

		self.g.serialize(destination=out_path, format="turtle")
		self.init_graph()

	def init_graph(self):
		if self.g is not None:
			self.g.remove((None, None, None))
			self.g.close()

		self.g = Graph()
		self.g.bind("gdsc", self.__gdsc_ns)
		self.g.bind("gdscc", self.__cellline_ns)
		self.g.bind("gdscd", self.__drug_ns)
		self.g.bind("m2r", self.__m2r_ns)
		self.g.bind("ontology", self.__ontology_ns)
		self.g.bind("dct", DCTERMS)
		self.g.bind("xsd", XSD)
		self.g.bind("sio", self.__sio_ns)
		self.g.bind("uo", self.__uo_ns)
		self.g.bind("tgo", self.__tgo_ns)
		self.g.bind("tgexpc", self.__tgexpc_ns)
		self.g.bind("bao", self.__bao_ns)
		self.g.bind("skos", self.__skos_ns)

	def create_gdsc_single_cell_turtle(self, is_debug=False):
		if is_debug:
			target_drug = ["Doxorubicin", "Etoposide", "Gemcitabine"]
			self.ic50 = self.ic50[self.ic50["Drug name"].isin(target_drug)]

		cosmic = str(self.ic50.loc[:, 'Cosmic sample Id'].values[0])
		id_name = self.get_cell_line_id(cosmic)
		if id_name == 'NO_ID':
			return self.g

		self.create_cell_line_turtle(cosmic, id_name)

		for index, item in self.ic50.iterrows():
			self.create_experiment_turtle(cosmic, id_name)
			self.create_drug_turtle(cosmic, id_name, item)

		if not is_debug:
			self.save_turtle(id_name)

		return self.g

	def create_cell_line_turtle(self, cosmic, cell_id):
		bot_id = self.get_bto_id(cell_id)
		if not bot_id == 'NO ID':
			self.g.add((self.__cellline_ns[cell_id], RDF["type"], self.__uo_ns[bot_id]))

		self.g.add((self.__cellline_ns[cell_id], RDF["type"], self.__m2r_ns["CellLine"]))
		self.g.add((self.__cellline_ns[cell_id], RDFS["label"], Literal("{}".format(cosmic))))
		self.g.add((self.__cellline_ns[cell_id], self.__ontology_ns['tcga_classification'], Literal(self.ic50['TCGA classification'].values[0])))
		self.g.add((self.__cellline_ns[cell_id], self.__ontology_ns['tissue'], Literal(self.ic50['Tissue'].values[0])))
		self.g.add((self.__cellline_ns[cell_id], self.__ontology_ns['tissue_subtype'], Literal(self.ic50['Tissue sub-type'].values[0])))

	def create_experiment_turtle(self, cosmic, cell_id):
		self.g.add((self.__cellline_ns[cell_id], self.__m2r_ns["drug"], BNode(cosmic + 'experiment')))
		self.g.add((BNode(cosmic + 'experiment'), RDF["type"], self.__tgo_ns["ExperimentalCondition"]))
		self.g.add((BNode(cosmic + 'experiment'), self.__tgo_ns['testType'], Literal('in vitro')))
		self.g.add((BNode(cosmic + 'experiment'), RDF["type"], self.__bao_ns["ExperimentalSpecification"]))

	def create_drug_turtle(self, cosmic, cellosaurus, item):
		self.g.add((BNode(cosmic + 'experiment'), self.__m2r_ns["drug"], self.__drug_ns[str(item['IC Result ID'])]))
		self.g.add((self.__drug_ns[str(item['IC Result ID'])], RDF["type"], self.__tgo_ns["ChemicalCompound"]))
		self.g.add((self.__drug_ns[str(item['IC Result ID'])], self.__skos_ns["altLabel"], Literal(cellosaurus.lower() + '_' + item['Drug name'].lower())))
		self.g.add((self.__drug_ns[str(item['IC Result ID'])], RDFS["label"], Literal(item['Drug name'])))
		self.g.add((self.__drug_ns[str(item['IC Result ID'])], self.__gdsc_ns['dataset_ver'], Literal(str(item['Dataset version']))))
		self.g.add((self.__drug_ns[str(item['IC Result ID'])], self.__gdsc_ns["ic50_results_id"], Literal(str(item['IC Result ID']))))
		self.g.add((self.__drug_ns[str(item['IC Result ID'])], self.__sio_ns["SIO_000216"], BNode((str(item['IC Result ID'])) + 'ic50')))
		self.g.add((BNode(str(item['IC Result ID']) + 'ic50'), self.__sio_ns['SIO_000300'], Literal(item['IC50'], datatype=XSD.decimal)))
		self.g.add((BNode(str(item['IC Result ID']) + 'ic50'), self.__sio_ns['SIO_000221'], self.__uo_ns["UO_0000064"]))
		self.g.add((BNode(str(item['IC Result ID']) + 'ic50'), RDF["type"], self.__uo_ns["MI_0641"]))
		self.g.add((self.__uo_ns["MI_0641"], RDFS['label'], Literal('IC50')))
		self.g.add((self.__uo_ns["UO_0000064"], RDFS['label'], Literal('μM')))
		self.g.add((self.__drug_ns[str(item['IC Result ID'])], self.__sio_ns["SIO_000216"], BNode(str(item['IC Result ID']) + 'auc')))
		self.g.add((BNode(str(item['IC Result ID']) + 'auc'), self.__sio_ns['SIO_000300'], Literal(item['AUC'], datatype=XSD.decimal)))
		self.g.add((BNode(str(item['IC Result ID']) + 'auc'), RDF["type"], self.__uo_ns["MCIT_C64774"]))
		self.g.add((self.__uo_ns["MCIT_C64774"], RDFS['label'], Literal('AUC')))

	def get_cell_line_id(self, cosmic):
		rtrn = self.id_map[self.id_map["GDSC"] == cosmic]["ID"]
		if len(rtrn) == 0:
			return 'NO_ID'
		return rtrn.values[0].replace(' ', '')

def process(row):
	gtr = GdscToRDF(row)
	gtr.create_gdsc_single_cell_turtle()
	del gtr

if __name__ == '__main__':
	with open(get_param(DbName.GDSC, 'in_gdsc_drug_file')[0], 'rb') as f:
		ic50 = pickle.load(f)

	jobs = int(get_param(DbName.COMMON, 'n_jobs')[0])
	cell_line_list = list(set(ic50['Cosmic sample Id']))
	print("GDSCデータからRDF/turtleファイルに変換")
	aprun = ParallelExecutor(n_jobs=jobs)
	aprun(total=len(cell_line_list))(delayed(process)(ic50[ic50['Cosmic sample Id'] == cosmic]) for cosmic in cell_line_list)
