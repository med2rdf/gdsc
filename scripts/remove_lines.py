import csv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--version', help='GDSC version. Current options are [1, 2]', type=int, required=True)
args = parser.parse_args()
version = args.version

csv_paths = [
    '../data/raw/PANCANCER_ANOVA.csv'.replace('.csv', '_GDSC' + str(version) + '.csv'),
    '../data/raw/PANCANCER_IC.csv'.replace('.csv', '_GDSC' + str(version) + '.csv')
]
for csv_path in csv_paths:
    with open(csv_path) as f:
        reader = csv.reader(f)
        header = next(reader)
        column_num = len(header)
        rows = []
        removed_rows = []
        for row in reader:
            if len(row) == column_num:
                rows.append(row)
            else:
                removed_rows.append(row)
        print('{}/{} lines removed'.format(len(removed_rows), len(rows) + len(removed_rows)))

    with open(csv_path, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)
