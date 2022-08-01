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
  parser.add_argument("-i", "--index", help="Index file", required=True)
  
  return parser.parse_args()

args = get_args()

def readRecords(fh1, fh2, fh3, fh4):
  records = [[fhx.readline().strip() for _ in range(4)] for fhx in [fh1, fh2, fh3, fh4]]
  if records[0][0] == "" or records[1][0] == "" or records[2][0] == "" or records[3][0] == "":
    return None
  return records

def addToCountDict(dict, value):
  if value not in dict:
    dict[value] = 1
  else:
    dict[value] += 1



indexes = set()
with open(args.index, 'r') as inFile:
  for line in inFile:
    line = line.strip().split("\t")
    indexes.add(line[-1])

indexBucketCounts = {}
with gzip.open(args.R1, "rt") as fh1, gzip.open(args.R2, "rt") as fh2, \
      gzip.open(args.R3, "rt") as fh3, gzip.open(args.R4, "rt") as fh4:
  while True:
    records = readRecords(fh1, fh2, fh3, fh4)
    if records == None:
      break #base condition, file has ended
        
    R1_record = records[0]
    I1_record = records[1]
    I2_record = records[2]
    R2_record = records[3]

    #these do not have trailing newlines
    I1Seq = I1_record[SEQ_INDEX]
    rcI2Seq = bioinfo.reverse_compliment(I2_record[SEQ_INDEX])
    newR1Header = f"{R1_record[HEADER_INDEX]} {I1Seq}-{rcI2Seq}"
    newR2Header = f"{R2_record[HEADER_INDEX]} {I1Seq}-{rcI2Seq}"

    if I1Seq not in indexes or rcI2Seq not in indexes:
      #at least one index is not in the provided list, unknown
      addToCountDict(indexBucketCounts, f"{I1Seq}-{rcI2Seq}")
      with open("unknown_R1.fastq", "w") as out_1, open("unknown_R2.fastq", "w") as out_2:
        out_1.write(f"{newR1Header}\n")
        out_2.write(f"{newR2Header}\n")
        for line in R1_record[1:]:
          out_1.write(f"{line}\n")
        for line in R2_record[1:]:
          out_2.write(f"{line}\n")
        continue
    
    if bioinfo.qual_score(I1_record[QUALITY_INDEX]) < QUAL_CUTOFF or bioinfo.qual_score(I2_record[QUALITY_INDEX]) < QUAL_CUTOFF:
      #at least one index is low quality, unknown
      addToCountDict(indexBucketCounts, f"{I1Seq}-{rcI2Seq}")
      with open("unknown_R1.fastq", "w") as out_1, open("unknown_R2.fastq", "w") as out_2:
        out_1.write(f"{newR1Header}\n")
        out_2.write(f"{newR2Header}\n")
        for line in R1_record[1:]:
          out_1.write(f"{line}\n")
        for line in R2_record[1:]:
          out_2.write(f"{line}\n")
        continue
    
    if I1Seq in indexes and rcI2Seq in indexes and I1Seq != rcI2Seq:
      #good indexes, hopped
      addToCountDict(indexBucketCounts, f"{I1Seq}-{rcI2Seq}")
      with open("hopped_R1.fastq", "w") as out_1, open(f"hopped_R2.fastq", "w") as out_2:
        out_1.write(f"{newR1Header}\n")
        out_2.write(f"{newR2Header}\n")
        for line in R1_record[1:]:
          out_1.write(f"{line}\n")
        for line in R2_record[1:]:
          out_2.write(f"{line}\n")
        continue

    if I1Seq in indexes and rcI2Seq in indexes and I1Seq == rcI2Seq:
      #good indexes, matched
      addToCountDict(indexBucketCounts, f"{I1Seq}-{rcI2Seq}")
      with open(f"{I1Seq}_R1.fastq", "w") as out_1, open(f"{I1Seq}_R2.fastq", "w") as out_2:
        out_1.write(f"{newR1Header}\n")
        out_2.write(f"{newR2Header}\n")
        for line in R1_record[1:]:
          out_1.write(f"{line}\n")
        for line in R2_record[1:]:
          out_2.write(f"{line}\n")
        continue
    raise Exception("Something not handled! PANIC!!!!!11!1!!1!1!1!!1!!1")

with open("stats.tsv", 'w') as outFile:
  indexCounts = sorted(indexBucketCounts.items(), key=lambda item: item[1])
  for item in indexCounts:
    outFile.write(f"{item[0]}\t{item[1]}\n")