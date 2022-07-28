# Demux Pseudocode
## Requirements
### Inputs
- 4 fastq files
- 1 text file with barcodes

### Outputs
- paired fastqs for real indexes (48 files)
- paired fastqs for reads with index hopping (2 files)
- paired fastqs for reads with incorrect or missing indexes (2 files)
- summary statistics (one tsv file) including 
    - number of reads with each actual index
    - number of reads with hopped indexes
    - number of reads with bad indexes

## Test input and output
### Inputs
Test input fastq files are in the `Test-input_FASTQ/` directory. Each file contains 3 records, they are based on the first 3 records in the full data.
The indexes of the first record have been modified to have the first index listed in the index file.
The second record will have the last and first indexes.
The third record will have a matched index not found in the index file.

Quality scores for all bases have been set to J for convenience.

### Outputs
Test output fastq files are in the `Test-output_FASTQ/` directory.
Naming is a best guess for now, and subject to change. The most likely changes would be using an identifier instead of an index for the matched files, and capitalization or shortening of the names of the not good output files.

The tsv file has not been made yet, I'm not sure if I can put anything other than fastqs in that directory.

## Pseudocode
### General strategy
- Read the index file to get the list of proper indexes.
- Open all input fastq files as text, not binary
- Read the input fastqs simultaneously, one record at a time
    - Reading will be done in the function readRecords
    - If readRecords returns None, break
    - Look at the indexes first to see which bucket the reads belong in. (I2/R3 must be reverse complimented) The bucket will be the index, "hopped", or "unknown"
        - If quality of either index is too low, reads belong in unknown
        - If either index contains N, reads belong in unknown
        - If I1 and RC(I2) are the same, and that value is in the list of proper indexes, reads belong in that index's file
        - If I1 and RC(I2) are not the same, but are both in the list of proper indexes, reads belong in the hopped file.
    - *Check the R1 and R2 quality to ensure that they BOTH pass quality filter? Or are they always output and the downstream user deals with low quality reads?*
    - Open R1 and R2 output files as specified by the bucket for appending
    - Append indexes to headers
    - *Modify R2 to show 2:N:0:1 instead of 4:N:0:1?*
    - Add to count dictionary of index pairs
- Write tsv in format `index1-RCindex2  count` (Sort dictionary by values)

### Functions
- `readRecords(fh1, fh2, fh3, fh4)`
    - Given file handles for R1, R2, R3, and R4, reads next 4 lines in each file and returns those lines as a 2D list or signals that a file has ended.
    - Inputs: file handles for each input fastq file
    - Output: list of lists containing [[header, seq, +, qual], for next record in all 4 files] OR None if any file has ended
    - Example output (EXTREMELY simplified): `[["@", "A", "+", "#"], ["@", "T", "+", "#"], ["@", "C", "+", "#"], ["@", "A", "+", "#"]]`
    - Return statement: `return [record1, record2, record3, record4]`  OR `return None`
- `reverseCompliment(seq)`
    - Returns the reverse compliment of an input sequence, N is retained as N.
    - How it could work:
        - reverse string
        - use translate on the reversed string to turn A into T, etc.
    - Input: String containing sequence to be reverse complimented
    - Output: String containing the reverse compliment.
    - Example input: "ATGCATGCN"
    - Example output: "NGCATGCAT"
    - Return statement: `return RC`