#!/usr/bin/env python
# Author: <YOU> <optional@email.address>

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "0.4"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNAbases = set('ATGCNatcgn')
RNAbases = set('AUGCNaucgn')

def reverse_compliment(seq: str) -> str:
  seq = seq[::-1]
  x = "GCAT"
  y = "CGTA"
  mytable = seq.maketrans(x, y)
  return seq.translate(mytable)

def convert_phred(letter: str) -> int:
  """Converts a single character into a phred score, only works with Phred+33 for now"""
  return ord(letter) - 33

def convert_binary_phred(qualAscii: int) -> int:
  """Converts a single binary character into a phred score, only Phred+33 for now"""
  return qualAscii - 33

def qual_score(phred_score: str) -> float:
  """Input: a string containing Phred+33 quality scores
  Returns the average score represented by that string"""
  rolling_sum = 0
  for letter in phred_score:
    rolling_sum += convert_phred(letter)
  return(rolling_sum / len(phred_score))

def validate_base_seq(seq,RNAflag=False):
  '''This function takes a string. Returns True if string is composed
  of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
  return set(seq)<=(RNAbases if RNAflag else DNAbases)

def gc_content(DNA):
  '''Input: string
  This function returns the GC content of an input string.
  It returns None if the string is not DNA'''
  DNA = DNA.upper()
#     if not validate_DNA_seq(DNA):
#        return None
  gc_count = DNA.count("G") + DNA.count("C")
  
  return(gc_count / len(DNA))

def oneline_fasta(fnameIN, fnameOUT):
  '''Inputs:  fnameIN: fasta input filename
              fnameOUT: fasta output filename (will be created or overwritten)
  This function is given the filename of a fasta file and the desired output filename.
  It writes the contents of the input file in 2-line fasta format in the output file.'''
  with open(fnameIN, 'r') as inFile, open(fnameOUT, 'w') as outFile:
    header = ""
    seq = ""
    for line in inFile:
      line = line.strip()
      if line.startswith(">"):
        if header != "": #Doesn't break on the first record
          #save the previous record before reading anything
          outFile.write(f"{header}\n{seq}\n")
          seq = "" #reset seq for the new record
        header = line #only get the gene chunk of the header, but drop the "gene:"
      else:
          seq += line
    else:
      outFile.write(f"{header}\n{seq}\n")


if __name__ == "__main__":
  # write tests for functions above
  assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
  assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
  assert validate_base_seq("Hi there!") == False, "Validate base seq fails to recognize nonDNA"
  assert validate_base_seq("Hi there!", True) == False, "Validate base seq fails to recognize nonDNA"
  print("Passed DNA and RNA tests")

  assert gc_content("GCGCGC") == 1
  assert gc_content("AATTATA") == 0
  assert gc_content("GCATGCAT") == 0.5
  print("correctly calculated GC content")
