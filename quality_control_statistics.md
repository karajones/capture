# Quality control and statistics

My primary goals when curating data for onward analyses:
1. Make sure SNPs come from a contiguous region of high coverage/depth (not just one small high coverage area or areas with a mix of high/low depth)
2. Remove loci with indels (which are difficult to call accurately) or high numbers of SNPs (see note below)
3. Remove potential duplicate loci/repeating regions

## Basic mapping statistics

Once the final bam has been output, all of the information about duplicates, unmapped reads, etc. has been removed, so if you want to compare the before and after for some basic read mapping stats, then this needs to be run on `.marked.bam` files.
> Note: You can get a lot more detailed output using `stats` rather than `flagstats` but since we’re just mapping back to short sequences rather than trying to align to a chromosome, I find `flagstats` works just fine.
```
for f in *.marked.bam; do samtools flagstat -O tsv $f > ../stats/${f%%.*}.flagstats.tsv; done
```

Example `flagstats.tsv` output with my notes below:
```
7901300	0	total (QC-passed reads + QC-failed reads)
7698926	0	primary
0	0	secondary
202374	0	supplementary
359148	0	duplicates
359148	0	primary duplicates
2480857	0	mapped
31.40%	N/A	mapped %
2278483	0	primary mapped
29.59%	N/A	primary mapped %
...
```

- primary: “best” alignment
- secondary: alternative alignment (e.g., second “best” alignment for the same primary read)
- supplementary: single read split and aligned to more than one site
- duplicates: PCR and optical duplicates
	- “duplicates” includes secondary/supplementary alignments while “primary duplicates” doesn’t
	- I had a lot of duplicates in this run because the library was overamplified, which is not ideal
- mapped %: total reads/mapped
- primary mapped %: primary/primary mapped

Here's a table with an example before/after:
| flag | marked bam | final bam |
| ----- | ----- | ----- |
| total (QC-passed reads + QC-failed reads)	| 68351828 | 2394757 |
| primary	| 68233758| 2394757 |
| secondary	| 0	| 0 |
| supplementary	| 118070 | 0 |
| duplicates	| 3094704	| 0 |
| primary duplicates	| 3094704	| 0 |
| mapped	| 6888788	| 2394757 |
| mapped %	| 10.08%	| 100.00% |
| primary mapped	| 6770718	| 2394757 |
| primary mapped %	| 9.92%	| 100.00% |

All the secondary, supplementary and duplicates were removed and only primary mapped reads remain.

## Find contiguous regions of high coverage/depth
