import csv
import numpy as np

p = 0.25

raw_bench = csv.reader(open('sample_code/submission_hour_chunk_means.csv'))

raw_ours = csv.reader(open('submissions/2012_04_29-10_10.csv'))

blended = csv.writer(open('submissions/blended25.csv','w'))

blended.writerow(raw_ours.next())

raw_bench.next()

for bRow,oRow in zip(list(raw_bench),list(raw_ours)):
	row = []
	row.extend(oRow[:5])
	row.extend([str((1-p)*float(ours)+p*float(bench)) for ours,bench in zip(oRow[5:],bRow[5:])])
	blended.writerow(row)
