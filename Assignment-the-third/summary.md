# Demux output summary

All fastqs are gzipped and stored in the `Assignment-the-third/out/` directory. There are also three summary tsvs.

**Important notes: I ran 1 base error correction on indexes, and labeled any index with average quality < 30 as unknown.**

## Samples

The count and percentage of reads that went into each file is stored in `out/stats_files.tsv`
| File basename | Record Count | Percentage of all records |
|---------------|--------------|---------------------------|
| TACCGGAT      | 71376106     | 19.64948315364762         |
| TCTTCGAC      | 40223838     | 11.07342038463195         |
| CTCTGGAT      | 32848576     | 9.043047833588924         |
| CTAGCTCA      | 16522520     | 4.548566692554029         |
| TGTTCCGT      | 15146782     | 4.169832937383457         |
| TCGAGAGT      | 10936729     | 3.010826511627145         |
| AGAGTCCA      | 10640614     | 2.9293075407821627        |
| TATGGCAC      | 10459887     | 2.8795543062486164        |
| TAGCCATG      | 10054938     | 2.768073882343361         |
| ATCATGCG      | 9481052      | 2.610085951632848         |
| AACAGCGA      | 8396171      | 2.311423666340731         |
| GTCCTAAG      | 8365105      | 2.3028713527184217        |
| AGGATAGC      | 8322630      | 2.2911781987524265        |
| GTAGCGTA      | 7639445      | 2.1031008028193288        |
| ACGATCAG      | 7611984      | 2.0955409275736505        |
| GCTACTCT      | 6753597      | 1.859231301831247         |
| ATCGTGGT      | 6544042      | 1.801541863824323         |
| GATCAAGG      | 6244848      | 1.7191752597583567        |
| CGATCGAT      | 5344118      | 1.4712088189863564        |
| CGGTAATC      | 4615765      | 1.270696899725747         |
| TCGGATTC      | 4264662      | 1.1740400089212089        |
| CACTTCAC      | 3927043      | 1.081095195528736         |
| TCGACAAG      | 3642158      | 1.0026677872273235        |
| GATCTTGC      | 3506814      | 0.9654082644404223        |
| hopped        | 613814       | 0.1689799083810072        |
| unknown       | 49763497     | 13.699640548730605        |

## Index swapping

Reads per index pair are summarized in `out/stats_indexes.tsv`. Matched index pairs and hopped index pairs are shown individually, but unknown pairs are grouped together. The contents of the output file will not be shown here, as there are too many entries.

The number of records that had hopped indexes (with index error correction of <= 1 base) was 613814, 0.1689799083810072% of the records

## Index error correcting

Error correcting is summarized in `out/stats_corrections.tsv`. The records shown for unknown are only records which failed due to quality scores and were corrected.
This avoids counting indexes which were too incorrect to be possibly salvagable, and could suggest that the quality cutoff could be lowered.

| File basename | Record count | Percentage of records with corrected index(es) |
|---------------|--------------|------------------------------------------------|
| unknown       | 9182632      | 53.48694325741607                              |
| TACCGGAT      | 2069033      | 12.051691788228183                             |
| TCTTCGAC      | 1074690      | 6.259848271096182                              |
| CTCTGGAT      | 685227       | 3.991306377893554                              |
| TGTTCCGT      | 359914       | 2.0964250440995182                             |
| CTAGCTCA      | 359625       | 2.0947416785240063                             |
| TCGAGAGT      | 278517       | 1.6223042560374576                             |
| TATGGCAC      | 264082       | 1.5382233491775505                             |
| AGAGTCCA      | 262248       | 1.5275406762865864                             |
| AGGATAGC      | 244573       | 1.4245874356389345                             |
| AACAGCGA      | 217980       | 1.2696886787199526                             |
| ATCATGCG      | 216437       | 1.2607010209932579                             |
| TAGCCATG      | 202680       | 1.1805693247222682                             |
| GTCCTAAG      | 200882       | 1.1700963444289452                             |
| GTAGCGTA      | 189244       | 1.1023073874469158                             |
| ATCGTGGT      | 186386       | 1.0856601251119233                             |
| ACGATCAG      | 170263       | 0.9917469653403765                             |
| GATCAAGG      | 158933       | 0.9257520450270584                             |
| GCTACTCT      | 142740       | 0.8314311496489862                             |
| CGATCGAT      | 118342       | 0.6893178163917635                             |
| CGGTAATC      | 117629       | 0.6851647379995839                             |
| TCGGATTC      | 101348       | 0.5903312607161654                             |
| hopped        | 96202        | 0.560356868842173                              |
| TCGACAAG      | 93617        | 0.5452997753726296                             |
| CACTTCAC      | 93403        | 0.5440532693755378                             |
| GATCTTGC      | 81361        | 0.47391109546441906                            |