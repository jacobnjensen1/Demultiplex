#!/usr/bin/env python

import bioinfo
import gzip
import argparse
import matplotlib.pyplot as plt

def get_args():
  parser = argparse.ArgumentParser(description="A program to visualize quality distribution")
  parser.add_argument("-l", "--length", help="length of reads in file", required=True, type=int)
  parser.add_argument("-f", "--file", help="input fastq file", required=True)
  parser.add_argument("-o", "--output_name", help="output file base name (eg. <output_name>.png)")
  return parser.parse_args()

args = get_args()

sumScores = [0 for _ in range(args.length)]
numRecords = 0
isBinary = None
with gzip.open(args.file, "rb") as inFile:
  for i, line in enumerate(inFile):
    if i % 4 == 3:
      if isBinary == None:
        if type(line) == bytes:
          isBinary = True
          print("IT'S BINARY")
        else:
          isBinary = False
          print("IT'S NOT BINARY")
      numRecords += 1
      if numRecords % 1000000 == 0:
        print(f"on record: {numRecords}")

      line = line.strip()
      #print(line)
      for j, char in enumerate(line):
        if isBinary:
          sumScores[j] += bioinfo.convert_binary_phred(char)
        else:
          sumScores[j] += bioinfo.convert_phred(char)
      
print(sumScores)
sumScores = [val / numRecords for val in sumScores]
print(sumScores)

plt.bar(x=range(args.length), height=sumScores, width=1)
plt.xlabel("position in read")
plt.ylabel("quality score")
plt.title("average quality score at each position")
plt.savefig(f"{args.output_name}.png")