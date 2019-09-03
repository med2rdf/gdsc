import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from rdflib import Literal, Namespace, Graph, BNode
from rdflib.namespace import RDF, RDFS, DCTERMS, XSD
import pickle
from scripts.omics_to_rdf import *
from joblib import delayed
import re
import math
import argparse

tag_list = ["Cell line name",
            "Drug name",
            "Drug Id",
            "IC50",
            "Cosmic sample Id",
            "TCGA classification",
            "Tissue",
            "Tissue sub-type",
            "AUC",
            "Dataset version",
            "Max conc",
            "target_pathway"]


class GdscToRDF(OmicsToRDF):
    def __init__(self, ic50_row, drugHash, version):
        super().__init__()
        self.version = version
        self.__result_path = get_param(DbName.GDSC, 'out_folder')[0]
        self.__cellline_ns = Namespace("https://www.cancerrxgene.org/translation/CellLine/")
        self.__drug_ns = Namespace("https://www.cancerrxgene.org/translation/Drug/")
        self.__cellosaurus_ns = Namespace("https://web.expasy.org/cgi-bin/cellosaurus/")
        self.__m2r_ns = Namespace("http://med2rdf/org/ontology/med2rdf#")
        self.__ontology_ns = Namespace("https://purl.jp/bio/10/gdsc/ontology#")
        self.__sio_ns = Namespace("http://semanticscience.org/resource/")
        self.__obo_ns = Namespace("http://purl.obolibrary.org/obo/")
        self.__tgo_ns = Namespace("http://purl.jp/bio/101/opentggates/ontology/")
        self.__bao_ns = Namespace("http://www.bioassayontology.org/bao#")
        self.__skos_ns = Namespace("http://www.w3.org/2004/02/skos/core#")

        self.__id_map = None
        self.__ic50 = ic50_row
        self.__anova = ic50_row
        self.drugHash = drugHash

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

            with open(get_param(DbName.GDSC, 'in_gdsc_anova_file')[0], 'rb') as f:
                anova = pickle.load(f)

            self.drugHash = {}
            for index, item in self.__ic50.iterrows():
                drug_id = 'gdsc' + str(self.version) + '_' + str(item['Drug Id'])
                drug_name = item['Drug name']
                self.drugHash[drug_id] = [drug_name, None]
            for drug_id, [drug_name, target_pathway] in self.drugHash.items():
                target = list(anova[anova['drug_name'] == drug_name]['target_pathway'].unique())
                if len(target) > 1:
                    target.remove('Unclassified')

                self.drugHash[drug_id][1] = target[0]
        return self.__ic50

    @ic50.setter
    def ic50(self, value):
        self.__ic50 = value

    def save_turtle(self, name):
        name = re.sub(r'\W', "_", name)
        out_path = self.__result_path + GDSC_PREFIX + name + EXT
        # print(out_path)
        idx = 2
        while(os.path.isfile(out_path)):
            out_path = self.__result_path + GDSC_PREFIX + name + str(idx) + EXT
            # print(out_path)
            idx += 1

        self.g.serialize(destination=out_path, format="turtle")
        self.init_graph()

    def init_graph(self):
        if self.g is not None:
            self.g.remove((None, None, None))
            self.g.close()

        self.g = Graph()
        self.g.bind("gdscc", self.__cellline_ns)
        self.g.bind("gdscd", self.__drug_ns)
        self.g.bind("m2r", self.__m2r_ns)
        self.g.bind("ontology", self.__ontology_ns)
        self.g.bind("dct", DCTERMS)
        self.g.bind("xsd", XSD)
        self.g.bind("sio", self.__sio_ns)
        self.g.bind("obo", self.__obo_ns)
        self.g.bind("tgo", self.__tgo_ns)
        self.g.bind("bao", self.__bao_ns)
        self.g.bind("skos", self.__skos_ns)

    def create_gdsc_single_cell_turtle(self, is_debug=False):
        if is_debug:
            target_drug = ["Doxorubicin", "SN-38", "Omipalisib"]
            self.ic50 = self.ic50[self.ic50["Drug name"].isin(target_drug)]

        cosmic = str(self.ic50.loc[:, 'Cosmic sample Id'].values[0])
        id_name = self.get_cell_line_id(cosmic)
        if id_name == 'NO_ID':
            return self.g

        self.create_cell_line_turtle(cosmic, id_name)
        for drug_id, [drug_name, target_pathway] in self.drugHash.items():
            self.create_drug_turtle(cosmic, drug_id, drug_name, target_pathway)
        for _, item in self.ic50.iterrows():
            self.create_experiment_turtle(cosmic, item)
        if not is_debug:
            self.save_turtle(id_name)

        return self.g

    def create_cell_line_turtle(self, cell_id, cell_name):
        bot_id = self.get_bto_id(cell_id)
        if not bot_id == 'NO ID':
            self.g.add((self.__cellline_ns[cell_id], RDF["type"], self.__obo_ns[bot_id]))

        self.g.add((self.__cellline_ns[cell_id], RDF["type"], self.__m2r_ns["CellLine"]))
        self.g.add((self.__cellline_ns[cell_id], RDFS["label"], Literal(cell_name)))
        self.g.add((self.__cellline_ns[cell_id], DCTERMS['identifier'], Literal(cell_id)))
        self.g.add((self.__cellline_ns[cell_id], self.__ontology_ns['tcga_classification'], Literal(self.ic50['TCGA classification'].values[0])))
        self.g.add((self.__cellline_ns[cell_id], self.__m2r_ns['site_primary'], Literal(self.ic50['Tissue'].values[0])))
        self.g.add((self.__cellline_ns[cell_id], self.__m2r_ns['site'], Literal(self.ic50['Tissue sub-type'].values[0])))

    def create_drug_turtle(self, cellosaurus, drug_id, drug_name, target_pathway):
        self.g.add((self.__drug_ns[drug_id], RDF["type"], self.__m2r_ns["Drug"]))
        self.g.add((self.__drug_ns[drug_id], RDFS["label"], Literal(drug_name)))
        # self.g.add((self.__drug_ns[drug_id], self.__drug_ns["chebi_id"], Literal(libchebipy.search(drug_name, True)[0].get_id()[len('CHEBI:'):])))
        if target_pathway:
            self.g.add((self.__drug_ns[drug_id], self.__ontology_ns['target_pathway'], Literal(target_pathway)))

    def create_experiment_turtle(self, cell_id, item):
        drug_id = 'gdsc' + str(self.version) + '_' + str(item['Drug Id'])
        ic50_value = math.exp(item['IC50'])
        self.g.add((self.__cellline_ns[cell_id], self.__m2r_ns["has_assay"], BNode(str(item['tmp_id']))))
        self.g.add((BNode(str(item['tmp_id'])), RDF["type"], self.__ontology_ns["Assay"]))
        self.g.add((BNode(str(item['tmp_id'])), self.__bao_ns['BAO_0002854'], self.__bao_ns['BAO_0020008']))
        self.g.add((BNode(str(item['tmp_id'])), self.__m2r_ns["sample"], self.__cellline_ns[cell_id]))
        self.g.add((BNode(str(item['tmp_id'])), self.__m2r_ns["drug"], self.__drug_ns[drug_id]))
        self.g.add((BNode(str(item['tmp_id'])), self.__ontology_ns['dataset_ver'], Literal(str(item['Dataset version']))))
        self.g.add((BNode(str(item['tmp_id'])), DCTERMS['identifier'], Literal(str(item['tmp_id']))))
        self.g.add((BNode(str(item['tmp_id'])), self.__sio_ns["SIO_000216"], BNode(str(item['tmp_id']) + 'res')))
        self.g.add((BNode(str(item['tmp_id']) + 'res'), self.__sio_ns['SIO_000300'], Literal(ic50_value, datatype=XSD.decimal)))
        self.g.add((BNode(str(item['tmp_id']) + 'res'), self.__sio_ns['SIO_000221'], self.__obo_ns["UO_0000064"]))
        self.g.add((BNode(str(item['tmp_id']) + 'res'), RDF["type"], self.__bao_ns["BAO_0000190"]))
        self.g.add((BNode(str(item['tmp_id'])), self.__sio_ns["SIO_000216"], BNode(str(item['tmp_id']) + 'auc')))
        self.g.add((BNode(str(item['tmp_id']) + 'auc'), self.__sio_ns['SIO_000300'], Literal(item['AUC'], datatype=XSD.decimal)))
        self.g.add((BNode(str(item['tmp_id']) + 'auc'), RDF["type"], self.__bao_ns["BAO_0002120"]))
        self.g.add((BNode(str(item['tmp_id'])), self.__sio_ns["SIO_000216"], BNode((str(item['tmp_id'])) + 'max_conc')))
        self.g.add((BNode(str(item['tmp_id']) + 'max_conc'), self.__sio_ns['SIO_000300'], Literal(item['Max conc'], datatype=XSD.decimal)))
        self.g.add((BNode(str(item['tmp_id']) + 'max_conc'), RDF["type"], self.__sio_ns["SIO_001088"]))

    def get_cell_line_id(self, cosmic):
        rtrn = self.id_map[self.id_map["GDSC"] == cosmic]["ID"]
        if len(rtrn) == 0:
            return 'NO_ID'
        return rtrn.values[0].replace(' ', '')


def process(row, drugHash, version):
    gtr = GdscToRDF(row, drugHash, version)
    gtr.create_gdsc_single_cell_turtle()
    del gtr


if __name__ == '__main__':
    version = args.version

    pkl_path = get_param(DbName.GDSC, 'in_gdsc_drug_file')[0].replace('.pkl', '_GDSC' + str(version) + '.pkl')
    with open(pkl_path, 'rb') as f:
        ic50 = pickle.load(f)

    pkl_path = get_param(DbName.GDSC, 'in_gdsc_anova_file')[0].replace('.pkl', '_GDSC' + str(version) + '.pkl')
    with open(pkl_path, 'rb') as f:
        anova = pickle.load(f)

    drugHash_pkl = '../data/drugHash.pkl'.replace('.pkl', '_GDSC' + str(version) + '.pkl')
    if os.path.exists(drugHash_pkl):
        with open(drugHash_pkl, 'rb') as f:
            drugHash = pickle.load(f)
    else:
        drugHash = {}
        for index, item in ic50.iterrows():
            drug_id = 'gdsc' + str(version) + '_' + str(item['Drug Id'])
            drug_name = item['Drug name']
            drugHash[drug_id] = [drug_name, None]
        for drug_id, [drug_name, target_pathway] in drugHash.items():
            target = list(anova[anova['drug_name'] == drug_name]['target_pathway'].unique())
            if len(target) == 0:
                continue
            if len(target) > 1:
                target.remove('Unclassified')
            drugHash[drug_id][1] = target[0]
        with open('../data/drugHash.pkl', 'wb') as f:
            pickle.dump(drugHash, f)

    jobs = int(get_param(DbName.COMMON, 'n_jobs')[0])
    cell_line_list = list(set(ic50['Cosmic sample Id']))
    print("GDSCデータからRDF/turtleファイルに変換")
    aprun = ParallelExecutor(n_jobs=jobs)
    aprun(total=len(cell_line_list))(delayed(process)(ic50[ic50['Cosmic sample Id'] == cosmic], drugHash, version) for cosmic in cell_line_list)
