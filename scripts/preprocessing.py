import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

import pickle
from scripts.omics_to_rdf import *


def make_id_map(text, db_name='GDSC'):
    id_txt, db_txt = None, None
    idx = 0
    df = pd.DataFrame([], columns=['ID', db_name])
    for txt in text:
        tmp = txt.strip('\n').split('   ')
        if tmp[0] == 'ID':
            id_txt = tmp[1]
        elif tmp[0] == 'DR':
            tmp[1] = tmp[1].split('; ')
            if tmp[1][0] == db_name:
                db_txt = tmp[1][1]

        if (tmp[0] == '//') & (db_txt is not None):
            idx += 1
            df.loc[str(idx)] = [id_txt, db_txt]
            id_txt, db_txt = None, None
    return df


def check_raw_file():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--version', help='GDSC version. Current options are [1, 2]', type=int, required=True)
    args = parser.parse_args()
    version = args.version

    csv_path = get_param(DbName.PRE, 'in_gdsc_drug_file')[0].replace('.csv', '_GDSC' + str(version) + '.csv')
    print(csv_path)
    pkl_path = get_param(DbName.GDSC, 'in_gdsc_drug_file')[0].replace('.pkl', '_GDSC' + str(version) + '.pkl')
    if not os.path.isfile(pkl_path):
        print("GDSC DRUGファイルをpickle化")
        df = pd.read_csv(csv_path, low_memory=False)
        df.drop_duplicates(['Drug name', 'Cell line name'], keep=False, inplace=True)
        df = df.fillna({'Tissue': 'unknown', 'Tissue sub-type': 'unknown', 'TCGA classification': 'unknown'})
        df['tmp_id'] = range(len(df))
        with open(pkl_path, 'wb') as f:
            pickle.dump(df, f)
        del df

    csv_path = get_param(DbName.PRE, 'in_gdsc_anova_file')[0].replace('.csv', '_GDSC' + str(version) + '.csv')
    pkl_path = get_param(DbName.GDSC, 'in_gdsc_anova_file')[0].replace('.pkl', '_GDSC' + str(version) + '.pkl')
    if not os.path.isfile(pkl_path):
        print("GDSC ANOVAファイルをpickle化")
        df = pd.read_csv(csv_path, low_memory=False)
        with open(pkl_path, 'wb') as f:
            pickle.dump(df, f)
        del df

    if not os.path.isfile(get_param(DbName.COMMON, 'in_cellosaurus_gdsc_file')[0]):
        print("cellosaurusファイルをgdscについてpickle化")
        with open(get_param(DbName.PRE, 'in_cellosaurus_file')[0]) as f:
            text = f.read().split('\n')
        df = make_id_map(text, 'GDSC')

        with open(get_param(DbName.COMMON, 'in_cellosaurus_gdsc_file')[0], 'wb') as f:
            pickle.dump(df, f)
        del df


check_raw_file()
