#!/usr/bin/env python

import bioinfo
import argparse
import gzip

HEADER_INDEX = 0
SEQ_INDEX = 1
QUALITY_INDEX = 3
QUAL_CUTOFF = 30

indexes = {"GTAGCGTA","CGATCGAT","GATCAAGG","AACAGCGA","TAGCCATG",\
  "CGGTAATC","CTCTGGAT","TACCGGAT","CTAGCTCA","CACTTCAC","GCTACTCT",\
  "ACGATCAG","TATGGCAC","TGTTCCGT","GTCCTAAG","TCGACAAG","TCTTCGAC",\
  "ATCATGCG","ATCGTGGT","TCGAGAGT","TCGGATTC","GATCTTGC","AGAGTCCA","AGGATAGC"}
#{} is both dictionary and set, these don't have :, so python knows it is a set

def get_args():
  parser = argparse.ArgumentParser(description="A program to summarize and visualize kmer distribution")
  parser.add_argument("-1", "--R1", help="gzipped R1 input fastq file", required=True)
  parser.add_argument("-2", "--R2", help="gzipped R2 input fastq file", required=True)
  parser.add_argument("-3", "--R3", help="gzipped R3 input fastq file", required=True)
  parser.add_argument("-4", "--R4", help="gzipped R4 input fastq file", required=True)
  parser.add_argument("-o", "--out", help="output directory, MUST ALREADY EXIST", required=True)
  
  return parser.parse_args()

args = get_args()

def readRecords(fh1, fh2, fh3, fh4):
  records = [[fhx.readline().strip() for _ in range(4)] for fhx in [fh1, fh2, fh3, fh4]]
  if records[0][0] == "" or records[1][0] == "" or records[2][0] == "" or records[3][0] == "":
    return None
  return records

def addToCountDict(dict, key):
  if key not in dict:
    dict[key] = 1
  else:
    dict[key] += 1

def writeRecords(baseName, R1Header, R2Header, R1Record, R4Record):
  with open(f"{args.out}/{baseName}_R1.fastq", "a") as out_1, \
      open(f"{args.out}/{baseName}_R2.fastq", "a") as out_2:
    out_1.write(f"{R1Header}\n")
    out_2.write(f"{R2Header}\n")
    for line in R1Record[1:]:
      out_1.write(f"{line}\n")
    for line in R4Record[1:]:
      out_2.write(f"{line}\n")

indexBucketCounts = {}
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

    #ONLY FOR TESTING
    #TODO: CHECK THIS OUT, IT'S FOR TESTING
    if recordCount == 5000000:
      print("Test is over")
      break
    #END TEST
    #REMOVE LATER

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

    if I1Seq not in indexes or rcI2Seq not in indexes:
      #at least one index is not in the provided list: unknown
      addToCountDict(indexBucketCounts, f"{I1Seq}-{rcI2Seq}")
      writeRecords("unknown", newR1Header, newR2Header, R1_record, R2_record)
      unknownCount += 1
      continue
    
    if bioinfo.qual_score(I1_record[QUALITY_INDEX]) < QUAL_CUTOFF or bioinfo.qual_score(I2_record[QUALITY_INDEX]) < QUAL_CUTOFF:
      #at least one index is low quality: unknown
      addToCountDict(indexBucketCounts, f"{I1Seq}-{rcI2Seq}")
      writeRecords("unknown", newR1Header, newR2Header, R1_record, R2_record)
      unknownCount += 1
      poorQualCount += 1
      continue
    
    if I1Seq in indexes and rcI2Seq in indexes and I1Seq != rcI2Seq:
      #good indexes but not same: hopped
      addToCountDict(indexBucketCounts, f"{I1Seq}-{rcI2Seq}")
      writeRecords("hopped", newR1Header, newR2Header, R1_record, R2_record)
      hoppedCount += 1
      continue

    if I1Seq in indexes and rcI2Seq in indexes and I1Seq == rcI2Seq:
      #good indexes and same: matched
      addToCountDict(indexBucketCounts, f"{I1Seq}-{rcI2Seq}")
      writeRecords(I1Seq, newR1Header, newR2Header, R1_record, R2_record)
      matchedCount += 1
      continue

    #If the code gets to this point, a set of records have not been assigned to a bucket
    #This should not be possible, so exit and panic if it happens
    raise Exception("Something not handled! PANIC!!!!!11!1!!1!1!1!!1!!1")

with open(f"{args.out}/counts.tsv", 'w') as outFile:
  #output counts of all index pairs
  indexCounts = sorted(indexBucketCounts.items(), key=lambda item: item[1])
  outFile.write("Index pair\tCount\n")
  for item in indexCounts:
    outFile.write(f"{item[0]}\t{item[1]}\n")

#TODO: maybe look at all possible good indexes?
with open(f"{args.out}/stats.tsv", 'w') as outFile:
  #output counts of records in file groups
  percentUnknown = (unknownCount / recordCount) * 100
  percentHopped = (hoppedCount / recordCount) * 100

  matchedCounts = {f"{id}-{id}": indexBucketCounts[f"{id}-{id}"] for id in indexes if f"{id}-{id}" in indexBucketCounts}
  samplePercentages = {id: ((matchedCounts[id] / recordCount) * 100) for id in matchedCounts.keys()}
  samplePercentages = sorted(samplePercentages.items(), key=lambda item: item[1])

  outFile.write("Sample indexes or grouping\tPercentage of all reads\n")
  for item in samplePercentages:
    outFile.write(f"{item[0]}\t{item[1]}\n")
  outFile.write(f"hopped\t{percentHopped}\n")
  outFile.write(f"unknown\t{percentUnknown}\n")

print(f"""Finished demultiplexing {recordCount} records
    {matchedCount} were good quality and matched
    {hoppedCount} had good quality indexes that hopped
    {unknownCount - poorQualCount} had a weird index
    {poorQualCount} had a poor quality index""")