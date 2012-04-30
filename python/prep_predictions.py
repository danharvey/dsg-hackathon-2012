#!/usr/bin/python
import csv
import random

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

days =['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday', 'Sunday']
def computeDayOfWeek(sRow, tRow):
  dayIndex = days.index(tRow[4]) + (int(sRow[2]) - int(tRow[2]) + int(tRow[5])) /24
  print dayIndex
  dayIndex = dayIndex % 7
  print "%s, %i .. %i, %i" % (tRow[4], int(tRow[2]), int(sRow[2]), dayIndex)
  return days[dayIndex]

raw_data = csv.reader(open('imputedDataWithFactorLoadings/imputedDataSetWithFactorLoadings.txt', 'r'), delimiter=' ')
factors = loadChunks(raw_data)
  
raw_data = csv.reader(open('raw_data/SubmissionZerosExceptNAs.csv', 'r'))
raw_data.next() # Skip header
submission = loadChunks(raw_data)

raw_data = csv.reader(open('raw_data/TrainingData.csv', 'r'))
raw_data.next() # Skip header
training = loadChunks(raw_data)

output = csv.writer(open('training/submission_training.csv', "w"))

for (chunkId, sChunk) in submission.iteritems():
  if chunkId in training:
    tChunk = training[chunkId]
    fChunk = factors[chunkId]
  
    # 12 factorised features to use in prediction.
    lastFactorRow = fChunk[-1]
    f = lastFactorRow[46:]
    tRow = tChunk[0]
  
    for sRow in sChunk:
      hour = sRow[3] 
      month = sRow[4]
      day = computeDayOfWeek(sRow, tRow)
      output.writerow([month, hour, day] + f)
  else:
    for sRow in sChunk:
      hour = sRow[3] 
      month = sRow[4]
      day = "Monday"
      output.writerow([month, hour, day] + ["0.1"]*24)
  
