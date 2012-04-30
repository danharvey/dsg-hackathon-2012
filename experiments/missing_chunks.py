#!/usr/bin/python
import csv


missing = 0
chunk = []
currentChunkId = 1
raw_data = csv.reader(open('../raw_data/TrainingData.csv', 'r'))
raw_data.next() # Skip header
for row in raw_data:
  # Process chunk if it's changed
  chunkId = int(row[1])
  if chunkId != currentChunkId:
    if (192 - len(chunk)) != 0:
      missing += 1
      
      missing_segments=[]
      previousRowId = 0
      for sample in chunk:
        rowId = int(sample[2])
        if rowId != previousRowId+1:
          missing_segments.append("%i: %i to %i" %(rowId - previousRowId - 2, previousRowId+1, rowId-1))
        previousRowId = rowId
      print missing_segments
    
    currentChunkId = chunkId
    chunk = []
  
  chunk.append(row)

print missing