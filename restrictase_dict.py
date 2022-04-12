import csv
ER_d = {}
with open('C:/Users/user/Desktop/endonuclease.csv') as er_csv:
    csv_read = csv.reader(er_csv, delimiter=',')
    line_count = 0
    for row in csv_read:
        ER_d[row[1]] = row[0]
print(ER_d)
