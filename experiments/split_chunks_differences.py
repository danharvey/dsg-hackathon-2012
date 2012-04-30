#!/usr/bin/python
import csv
import random

ranges = [1, 2, 3, 4, 5, 10, 17, 24, 48, 72]

outputs = {}
for r in ranges:
  outputs["%i_train" % r] = csv.writer(open('training_diff/training_%i_train.csv' % r, "w"))
  outputs["%i_test" % r] = csv.writer(open('training_diff/training_%i_test.csv' % r, "w"))
  outputs["%i_ubertrain" % r] = csv.writer(open('training_diff/training_%i_ubertrain.csv' % r, "w"))

def processChunk(fChunk, tChunk):
  hours = {}
  for row in fChunk:
    hours[int(row[2])] = row
  
  tHours = {}
  for row in tChunk:
    tHours[int(row[2])] = row
  
  # Spit out rows for 1, 2, 3, 4, 5, 10, 17, 24, 48, and 72 hour after.
  for index,r1 in enumerate(fChunk[1:]):
    t1 = int(r1[2])
    for r in ranges:
      t2 = t1 + r
      if t2 in hours:
        r2 = hours[t2]
        x = r2[3:46] + r1[46:]
	x_prev = fChunk[index][46:]
        if (random.random() > 1.0/3.0):
          y = r2[46:]
          outputs["%i_train" % r].writerow(x + x_prev + y)
        else:
          y = tHours[t2][56:]
          outputs["%i_test" % r].writerow(x + x_prev + y)
        
        y = r2[46:]
        outputs["%i_ubertrain" % r].writerow(x + x_prev + y)

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

raw_data = csv.reader(open('imputedDataWithFactorLoadings/imputedDataSetWithFactorLoadings.txt', 'r'), delimiter=' ')
factors = loadChunks(raw_data)

raw_data = csv.reader(open('raw_data/TrainingData.csv', 'r'))
raw_data.next() # Skip header
training = loadChunks(raw_data)

for (chunkId, fChunk) in factors.iteritems():
  tChunk = training[chunkId]
  
  processChunk(fChunk, tChunk)
