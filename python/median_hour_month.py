#!/usr/bin/python
import numpy as np
import csv

def safefloat(s):
        try:
                return float(s)
        except: 
                return np.nan

def safemedian(m):
	n,d = m.shape
	result = np.zeros(d)
	num = np.zeros(d)
	error = np.zeros(d)
	for i in range(0,d):
		col = np.array(filter(lambda(x):~np.isnan(x),m[:,i]))
		num[i] = len(col)
		if num[i]!=0:
			result[i] = np.median(col)
			error[i] = np.mean(col-np.median(col))
		else:
			result[i] = np.nan
			error[i] = np.nan
	return result,num,error	

contents = csv.reader(open('raw_data/TrainingData.csv'))

grouped = dict()

for row in contents:
	if row[3]=='10' or row[3]=='9' or row[3]=='11':
        	try:
                	grouped[int(row[5])].append([safefloat(a) for a in row[95-39:95]])
        	except:
                	grouped[int(row[5])]=[[safefloat(a) for a in row[95-39:95]]]

med = dict()
num = dict()
err = dict()

for hour in grouped:
        med[hour],num[hour],err[hour] = safemedian(np.array(grouped[hour]))

hours = med.keys()
hours.sort()
for hour in hours:
	print '10,'+str(hour)+','+','.join([str(a) for a in med[hour]])
