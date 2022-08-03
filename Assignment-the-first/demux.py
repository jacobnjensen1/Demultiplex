#!/usr/bin/env python

import bioinfo
import argparse
import gzip

HEADER_INDEX = 0
SEQ_INDEX = 1
QUALITY_INDEX = 3
QUAL_CUTOFF = 30


def get_args():
  parser = argparse.ArgumentParser(description="A program to summarize and visualize kmer distribution")
  parser.add_argument("-1", "--R1", help="gzipped R1 input fastq file", required=True)
  parser.add_argument("-2", "--R2", help="gzipped R2 input fastq file", required=True)
  parser.add_argument("-3", "--R3", help="gzipped R3 input fastq file", required=True)
  parser.add_argument("-4", "--R4", help="gzipped R4 input fastq file", required=True)
  parser.add_argument("-i", "--indexFile", help="file containing used indexes", required=True)
  parser.add_argument("-o", "--out", help="output directory, MUST ALREADY EXIST", required=True)
  
  return parser.parse_args()

args = get_args()

# indexes_set = {"GTAGCGTA","CGATCGAT","GATCAAGG","AACAGCGA","TAGCCATG",\
#   "CGGTAATC","CTCTGGAT","TACCGGAT","CTAGCTCA","CACTTCAC","GCTACTCT",\
#   "ACGATCAG","TATGGCAC","TGTTCCGT","GTCCTAAG","TCGACAAG","TCTTCGAC",\
#   "ATCATGCG","ATCGTGGT","TCGAGAGT","TCGGATTC","GATCTTGC","AGAGTCCA","AGGATAGC"}
# {} is both dictionary and set, these don't have :, so python knows it is a set

#don't hardcode the indexes
indexes = set()
with open(args.indexFile) as inFile:
  inFile.readline()
  for line in inFile:
    line = line.strip().split()
    indexes.add(line[-1])

#Key= index, "hopped", or "unknown": value = tuple (file handle for R1, file handle for R2)
files = {index: (open(f"{args.out}/{index}_R1.fastq", 'w'), open(f"{args.out}/{index}_R2.fastq", 'w')) for index in indexes}
files["hopped"] = (open(f"{args.out}/hopped_R1.fastq", 'w'), open(f"{args.out}/hopped_R2.fastq", 'w'))
files["unknown"] = (open(f"{args.out}/unknown_R1.fastq", 'w'), open(f"{args.out}/unknown_R2.fastq", 'w'))

def readRecords(fh1, fh2, fh3, fh4):
  """Reads through each file to grab one record from each, returns the records in a 2D list[[R1Header, R1Seq, ...],...]"""
  records = [[fhx.readline().strip() for _ in range(4)] for fhx in [fh1, fh2, fh3, fh4]]
  if records[0][0] == "" or records[1][0] == "" or records[2][0] == "" or records[3][0] == "":
    return None
  return records

def addToCountDict(dict, key):
  """Given a key and a dictionary, increase the count for that key"""
  if key not in dict:
    dict[key] = 1
  else:
    dict[key] += 1

def writeRecords(baseName, R1Header, R2Header, R1Record, R4Record):
  """Uses basename to find output files from files dictionary, basename can be index1, 'hopped', or 'unknown'
  requires modified headers with indexes and full records. Writes with trailing newlines."""
  out_1,out_2 = files[baseName]
  out_1.write(f"{R1Header}\n")
  out_2.write(f"{R2Header}\n")
  for line in R1Record[1:]:
    out_1.write(f"{line}\n")
  for line in R4Record[1:]:
    out_2.write(f"{line}\n")

def errorCorrectIndex(qIndex):
  """Checks if an index is one character off from a proper index. 
  If there is a proper index that close, it returns that index, otherwise (or if there are more options),
  it returns the input index."""
  possibleIndexes = []
  for index in indexes:
    runSum = 0
    for i in range(len(index)):
      if qIndex[i] != index[i]:
        runSum += 1
      if runSum > 1:
        break
    if runSum <= 1:
      possibleIndexes.append(index)

  #if there's only one possible proper index: return that
  if len(possibleIndexes) == 1:
    return possibleIndexes[0]
  #if there is not only one possible index (0 or 2+): give up and return the input
  if len(possibleIndexes) != 1:
    return qIndex

indexBucketCounts = {}
correctionOutcomeCounts = {}
recordCount = 0
matchedCount = 0
hoppedCount = 0
unknownCount = 0
poorQualCount = 0

with gzip.open(args.R1, "rt") as fh1, gzip.open(args.R2, "rt") as fh2, \
      gzip.open(args.R3, "rt") as fh3, gzip.open(args.R4, "rt") as fh4:
  while True:
    records = readRecords(fh1, fh2, fh3, fh4)
    if records == None:
      print("Done reading")
      break #At least one file has ended

    if recordCount % 2000000 == 0:
      print(f"On Record {recordCount}")
    recordCount += 1

    R1_record = records[0]
    I1_record = records[1]
    I2_record = records[2]
    R2_record = records[3]

    #these do not have trailing newlines
    I1Seq = I1_record[SEQ_INDEX]
    rcI2Seq = bioinfo.reverse_compliment(I2_record[SEQ_INDEX])
    newR1Header = f"{R1_record[HEADER_INDEX]} {I1Seq}-{rcI2Seq}"
    #change read id from 4 to 2 for R4/R2
    newR2HeaderPortions = R2_record[HEADER_INDEX].split(" ")
    newR2HeaderDir = f"2{newR2HeaderPortions[1][1:]}"
    newR2Header = f"{newR2HeaderPortions[0]} {newR2HeaderDir} {I1Seq}-{rcI2Seq}"

    
    #Lesie's challenge for me:
    #error correct indexes, they can be off by one.
    #I think I want to do that here, before bucketing

    isCorrected = False
    if I1Seq not in indexes:
      I1Seq = errorCorrectIndex(I1Seq)
      isCorrected = True
    if rcI2Seq not in indexes:
      rcI2Seq = errorCorrectIndex(rcI2Seq)
      isCorrected = True

    #possible index correction summaries:
    #option 1
    #I1Q  I1C I2Q I2C Destination Count
    #AAA  AAT AAG AAT AAT 1
    #whatever whatever  whatever  whatever  unknown 10000 

    #option 2
    #Destination after correction  count percentage
    #AAT  12  x%
    #unkown 100000  a lot%

    #possible correction storage:
    #for option 1: {bucket: {(I1Q, I1C, I2Q, I2C), count}} #messy

    #for option 2: {bucket: count}  calculate percentage at end

    # #ONLY FOR TESTING		
    # #TODO: CHECK THIS OUT, IT'S FOR TESTING		
    # if recordCount == 5000000:		
    #   print("Test is over")		
    #   break		
    # #END TEST		
    # #REMOVE LATER

    if I1Seq not in indexes or rcI2Seq not in indexes:
      #at least one index is not in the provided list: unknown
      addToCountDict(indexBucketCounts, f"unknown")
      writeRecords("unknown", newR1Header, newR2Header, R1_record, R2_record)
      unknownCount += 1
      if isCorrected:
        addToCountDict(correctionOutcomeCounts, "unknown_no_correction")
      continue
    
    if bioinfo.qual_score(I1_record[QUALITY_INDEX]) < QUAL_CUTOFF or bioinfo.qual_score(I2_record[QUALITY_INDEX]) < QUAL_CUTOFF:
      #at least one index is low quality: unknown
      addToCountDict(indexBucketCounts, f"unknown")
      writeRecords("unknown", newR1Header, newR2Header, R1_record, R2_record)
      unknownCount += 1
      poorQualCount += 1
      if isCorrected:
        addToCountDict(correctionOutcomeCounts, "unknown_low_qual")
      continue
    
    if I1Seq != rcI2Seq:
      #good indexes but not same: hopped
      addToCountDict(indexBucketCounts, f"{I1Seq}-{rcI2Seq}")
      writeRecords("hopped", newR1Header, newR2Header, R1_record, R2_record)
      hoppedCount += 1
      if isCorrected:
        addToCountDict(correctionOutcomeCounts, "hopped")
      continue

    if I1Seq == rcI2Seq:
      #good indexes and same: matched
      addToCountDict(indexBucketCounts, f"{I1Seq}-{rcI2Seq}")
      writeRecords(I1Seq, newR1Header, newR2Header, R1_record, R2_record)
      matchedCount += 1
      if isCorrected:
        addToCountDict(correctionOutcomeCounts, I1Seq)
      continue

    #If the code gets to this point, a set of records have not been assigned to a bucket
    #This should not be possible, so exit and panic if it happens
    raise Exception("Something not handled! PANIC!!!!!11!1!!1!1!1!!1!!1")

with open(f"{args.out}/counts.tsv", 'w') as outFile:
  #output counts of all index pairs
  indexCounts = sorted(indexBucketCounts.items(), key=lambda item: item[1], reverse=True)
  outFile.write('Index pair or "unknown"\tRecord count\n')
  for item in indexCounts:
    outFile.write(f"{item[0]}\t{item[1]}\n")

with open(f"{args.out}/stats.tsv", 'w') as outFile:
  #output counts of records in file groups
  percentUnknown = (unknownCount / recordCount) * 100
  percentHopped = (hoppedCount / recordCount) * 100

  matchedCounts = {f"{id}": indexBucketCounts[f"{id}-{id}"] for id in indexes if f"{id}-{id}" in indexBucketCounts}
  samplePercentages = {id: ((matchedCounts[id] / recordCount) * 100) for id in matchedCounts.keys()}
  samplePercentages = sorted(samplePercentages.items(), key=lambda item: item[1], reverse=True)

  outFile.write("File basename\tPercentage of all reads\n")
  for item in samplePercentages:
    outFile.write(f"{item[0]}\t{item[1]}\n")
  outFile.write(f"hopped\t{percentHopped}\n")
  outFile.write(f"unknown\t{percentUnknown}\n")

with open(f"{args.out}/errorCorrectionOutcomes.tsv", 'w') as outFile:
  sumOutcomes = sum(correctionOutcomeCounts.values())
  outFile.write(f"File basename (and explanation)\tRecord count\tPercentage of records with corrected index(es)\n")
  for dest, count in sorted(correctionOutcomeCounts.items(), key=lambda item: item[1], reverse=True):
    outFile.write(f"{dest}\t{count}\t{(count/sumOutcomes) * 100}\n")

for fileTuple in files.values():
  fileTuple[0].close()
  fileTuple[1].close()

print(f"""Finished demultiplexing {recordCount} records
    {matchedCount} were good quality and matched
    {hoppedCount} had good quality indexes that hopped
    {unknownCount - poorQualCount} had a weird index
    {poorQualCount} had a poor quality index""")