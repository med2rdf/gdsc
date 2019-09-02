import csv

with open('data/raw/PANCANCER_ANOVA.csv') as f:
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

with open('data/raw/PANCANCER_ANOVA.csv', 'w') as f:
    writer = csv.writer(f)
    writer.writerow(header)
    writer.writerows(rows)
