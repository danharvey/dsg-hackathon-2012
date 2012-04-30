#!/usr/bin/python
import csv
import random
from datetime import datetime

def loadChunks(raw_data):
  chunks = {}

  currentChunkId = 1
  chunk = []

  for row in raw_data:
    chunkId = int(row[1])

    # Process chunk if it's changed
    if chunkId != currentChunkId:
      chunks[currentChunkId] = chunk

      chunk = []
      currentChunkId = chunkId

    chunk.append(row)

  chunks[currentChunkId] = chunk

  return chunks

raw_data = csv.reader(open('raw_data/TrainingData.csv', 'r'))
raw_data.next() # Skip header
training = loadChunks(raw_data)

dayPreds = {}
raw_data = csv.reader(open('missing_trainingchunk_predictions/median_3month.csv', 'r'))
for row in raw_data:
  dayPreds["%s_%s" % (row[0], row[1])] = row[2:]

sub = csv.reader(open('raw_data/SubmissionZerosExceptNAs.csv', 'r'))

output = csv.writer(open('submissions/%s.csv' % (datetime.now().strftime("%Y_%m_%d-%H_%M")), "w"))
output.writerow(sub.next())

predictions = csv.reader(open('training/finalResult.txt', 'r'), delimiter=' ')

t=0
for row in sub:
  pred = predictions.next()
  
  
  chunkId = int(row[1])
  # If not got a past chunk we can use...
  if chunkId not in training:
    hour = row[3] 
    month = row[4]
    pred = dayPreds["%s_%s" % (month, hour)]
  
  # Set -1,000,000 predictions...
  new_pred = []
  for index, sItem in enumerate(row[5:]):
    if sItem == "-1e+06":
      new_pred.append(sItem)
    else:
      new_pred.append(pred[index])
  
  
  for (s, p) in zip(row[5:], new_pred):
    diff = abs(float(s)-float(p))
    t += diff
    if diff > 100:
      print len(pred)
      print diff
      print row[5:]
      print new_pred
  
  output.writerow(row[0:5] + new_pred)

print t/(2100*39*1.0)
