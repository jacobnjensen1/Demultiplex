# Assignment the third notebook

I was challenged by Leslie to include error correction of index reads. If an input index read is within one base of a correct index, it gets changed to that correct index.

I wrote a function to check for that by comparing input indexes to all real indexes character by character.

Jason pointed out that it would be much more efficient to pre-calculate a dictionary of off-by-one indexes and check for those in the inputs.
I didn't do that because what I wrote was fast enough and it already worked.
If I really needed maximum speed, I would have implemented it. For instance, if we were demultiplexing over 1 billion reads, it would certainly make sense to do it.

I made demux.py in the `Assignment-the-first/` directory for some reason, so some slurm outputs are there.
I did a number of 5,000,000 record tests to make sure that my script was working and to clarify my outputs.
I copied the last slurm output into `Assignment-the-third/`.

The interesting parts of the last output:
```
Finished demultiplexing 363246735 records
    312869424 were good quality and matched
    613814 had good quality indexes that hopped
    13615974 had a weird index
    36147523 had a poor quality index
Command being timed: "python -u ./demux.py -o out -1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -3 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -4 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -i /projects/bgmp/shared/2017_sequencing/indexes.txt"
	User time (seconds): 4736.30
	System time (seconds): 48.25
	Percent of CPU this job got: 94%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:24:49
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 241816
	Average resident set size (kbytes): 0
	Exit status: 0
```

As can be seen, my script ran in 4736.20 seconds, or 1 hour and 19 minutes. I consider this fast enough, but I would certainly hope that mainstream tools have much shorter runtimes.

To keep the output directory small, I gzipped the fastq files after writing them out.
I decided not to gzip while writing because of suggestions from a stack overflow discussion that I unfortunately cannot find anymore.
The basic issue is that a lot of text is required to get a good compression ratio, and writing 4 lines out of millions at a time does not count as a lot of text.
(I think this might be due to including a large number of encoding trees in the output, but I could be wrong on that)
The person asking the question said that when they attempted it, their compressed file was about the same size as their uncompressed file.