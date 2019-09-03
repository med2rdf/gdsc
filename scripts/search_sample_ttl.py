from rdflib import Graph
import argparse
from glob import glob


def search_ttl(file_path):
    drugs = []
    with open(file_path) as f:
        while True:
            line = f.readline().strip()
            if line.startswith('gdscd'):
                drugs.append(line.split(' ')[0].split(':')[-1])
            if line.startswith('m2r:drug'):
                drug = line.split(' ')[1].split(':')[1]
                if drug in drugs:
                    print(file_path, drug)
                    return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--version', help='GDSC version. Current options are [1, 2]', type=int, required=True)
    args = parser.parse_args()
    version = args.version
    ttl_files = glob('../output/result_gdsc_{}*{}'.format(version, '.ttl'))

    for ttl_file in ttl_files:
        if search_ttl(ttl_file):
            break
