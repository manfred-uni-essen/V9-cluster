# V9-cluster
R-Script for modifying DNA-sequence-abundance tables: clustering of related sequences (e.g. SSU-ITS1) according to 100% identical sub-sequences; e.g. 150 BP of V9-SSU;

in our studies we use primers that include the V9 subregion of the SSU-ribosomal DNA plus the directly following ITS1-region;
whilest the V9 is largely conserved, the ITS1 region is more variable - 
this is the reason for V9 clustering - regarding identical V9 regions as OTU´s (operational taxonomic units) and variation in ITS1 as variation of these OTU´s - like species and subspecies in classical taxonomy.

What does the script :

In a first step the sequences are ordered according to their V9 subsequence, and the corresponding ordered abundance table is stored.

In a second part, the full sequences (V9 plus ITS1) of this ordered abundance table are transformed and stored in the fasta Format.
The fasta-formatted file can be read by several alignment programs allowing for manual assessment of the sequence data within the ITS1: the variation within the ITS1 is clearly visible by alignment.

In the third and last step, OTU´s with identical V9 are clustered, i.e. the abundances are summed up by aggregation. Of course, there is the need to choose one representative sequence for each V9 group; normally, this is the most frequent sequence - which is reasonable in most cases (default). Alternatively, the representative with the longest V9+ITS1 chain, i.e. sequence length can be selected. 
In any case, the resulting table is stored to a file - it is of course always smaller (less OTUs) than the original one - usually by a factor of 2. Within the result table, there is some additional information within the first columns like: 
V9 group: 1....n, plus the number of variants within each group.
