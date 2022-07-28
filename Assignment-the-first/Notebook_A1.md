# Assignment the first notebook

To view gzipped files: `zcat $file`
R1, R4 have long reads, they must be bio data.
R2, R3 have short reads, they must be indexes.
Headers are aligned across the 4 files. This will be helpful for parsing.

`zcat 1294_S1_L008_R1_001.fastq.gz | head -2 | tail -1 | wc -L` gives the length of the first sequence in R1: 101. That command ignores newlines.

Looking at the quality scores, there are many capital letters and #s. These are only seen in Phred+33.

Looking at the indexes themselves, the first many indexes in R2 have a leading N. If that N is removed, the remaining chunk matches the last 7 nucleotides in an actual index.
R3 indexes can also start with N, but the last 7 nucleotides are not found in actual indexes. If the non-N R3 indexes are reverse complimented, they can be found in the actual indexes.

`zcat 1294_S1_L008_R3_001.fastq.gz | sed -n 2~4p | grep -c N` counts the number of indexes in R3 that contain N, it does not multicount indexes with multiple Ns (tested with small test file that no longer exists).