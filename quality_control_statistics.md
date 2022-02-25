# Quality control and statistics

The ultimate goal here is to curate a list of loci and positions that will produce high quality sets of data for onward analyses. Quality control includes:
1. Making sure SNPs come from a contiguous region of high coverage/depth (not just one small high coverage area or areas with a mix of high/low depth)
2. Removing loci with indels (which are difficult to call accurately) or high numbers of SNPs (see note below)
3. Removing potential duplicate loci/repeating regions

In the first section, I'll go over looking at data quality and in the second section I'll cover how to implement quality controls.

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
### What do these flags mean?
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

All the secondary, supplementary and duplicates were removed and only primary mapped reads remain, which is good!

# Coverage statistics

>[Bedtools](https://bedtools.readthedocs.io/en/latest/) is required for most of these analyses.

Bedtools `genomecov` output statistics on the depth of reads at each individual site on the reference locus. I've combined it with `merge` to look at the coverage of reads mapped across each locus. Loci with no mapped reads are removed.
- `genomecov -bg`: output bed graph format with no zero values included
- `merge -d 500`: regions must be within 500 bp of each other (a number large enough to cover everything on one locus)
- `-c 4 -o min,max,mean,median`: calculate min, max, mean, and median for depth (column 4 on bedgraph output)
- the `echo` command is just used to add a header so the file is easier to read by humans

```
bedtools genomecov -bg -ibam DWR12.final.bam | bedtools merge -d 500 -i - -c 4 -o min,max,mean,median > DWR12.stats.bed
echo -e "locus\tstart\tend\tmin\tmax\tmean\tmedian\n$(cat DWR12.stats.bed)" > DWR12.stats.bed
```
Example output: (1) locus name, (2) start of coverage on the locus (if there is zero coverage on a locus, that region is omitted), (3) end of coverage on the locus, (4) minimum depth of coverage (on at least one position on the locus), (5) maximum depth of coverage, (6) mean depth of coverage across all positions on the locus, (7) median depth of coverage across the locus.
``` 
locus	start	end	min	max	mean	median
D01-LOCUS392-290	30	285	2	12	6.363636364	6
D01-LOCUS478-290	0	292	2	88	45.69230769	36
D01-LOCUS491-290	102	123	2	2	2	2
D01-LOCUS563-290	124	143	2	2	2	2
D01-LOCUS1024-290	29	176	2	4	2.5	2
D01-LOCUS1179-290	51	70	2	2	2	2
D01-LOCUS1313-290	0	298	2	14	7.033333333	6.5
D01-LOCUS1458-290	8	72	2	4	3	3
D01-LOCUS1675-290	0	143	4	188	122.7209302	149
...
```

## Create a histogram of coverage
Output the coverage at each depth across *all* loci (rather than per locus, as above). This can be used to see what the distribution of depths are across all loci and create a cumulative summation of depth to identify a cut-off for depth. The main `genomecov` output includes two sections: (1) information about read coverage for each region of each reference locus that has mapped reads and (2) depth across the entire "genome" (i.e., all reference loci combined). I only want the section for the entire genome , so I `grep` the output for that section.

```
bedtools genomecov -ibam DWR12.final.bam | grep ^genome > DWR12.coverage.hist.txt 
```

Or make all the files in one go:
```
for f in *.final.bam; do bedtools genomecov -ibam $f | grep ^genome > ${f%%.*}.coverage.hist.txt; done
```

Example output: (1) ignore this column, (2) depth, (3) number of bases at this depth, (4) count of all bases across all loci (will always be the same), (5) fraction of bases at this depth
```
genome	0	1424666	2100808	0.678151
genome	1	5283	2100808	0.00251475
genome	2	247435	2100808	0.117781
genome	3	3838	2100808	0.00182692
genome	4	82906	2100808	0.0394639
genome	5	3312	2100808	0.00157654
genome	6	42295	2100808	0.0201327
genome	7	2902	2100808	0.00138137
genome	8	28932	2100808	0.0137718
genome	9	2801	2100808	0.0013333
```
So, in the above example there are 1,424,666 bases with zero depth and 247,435 with a depth of two. What does that look like when graphed cumulatively?

>Note: This graph only shows bases that have at least one read mapped to them. I’ve set the x-axis of the graph a max depth of 200 but the depth goes up to 28,342!


## Find contiguous regions of high coverage/depth

### Output loci/regions with high depth
I use `bedtools genomecov` and `merge` functions to combine any continuous sites on a locus with at least 20 bp depth into a single region. The output is `bed` format.
1. `genomecov -bga`: calculate coverage per site (bp); exclude zero values
2. `sort`: data by locus (`k1,1`) and start position (`-k2,2n`)
3. `awk`: find positions with at least 20 bp depth
4. `merge`: merge all overlapping and abutting positions into a single contiguous feature

>I’m just running one bam file at a time here.

```
bedtools genomecov -bga -ibam H1.final.bam | sort -k1,1 -k2,2n | awk '$4 >= 20' | bedtools merge -i - > H1.clean.bed
```
Output columns: (1) locus name, (2) start position, (3) end position
```
D01-LOCUS100763-290	163	274
D01-LOCUS101326-290	43	130
D01-LOCUS101329-290	37	176
D01-LOCUS101329-290	193	287
D01-LOCUS10164-290	89	143
D01-LOCUS10164-290	207	227
D01-LOCUS10164-290	263	298
D01-LOCUS101812-290	192	298
D01-LOCUS1024-290	4	68
D01-LOCUS1024-290	77	86
...
```

### Print list of loci at least 100 contiguous bp long
Uses the `bed` file produced above to create a file which lists *one* region for each locus where all bp in a region are at least 100 bp long and 20 bp depth (or whatever was specified in the `bed` file). This command can be piped with the previous command to output this file without the intermediate `bed` file. Leave out the `,$3-$2` from the `print` command to output a file with just the positions. Leave out `-u` from `sort` command to output *all* regions.

```
cat H1.clean.bed | awk -F"\t" '$3-$2>=100{print $1,$2,$3,$3-$2}' | sort -u -k1,1 > H1.lengths
```

Output columns: (1) locus name, (2) start position, (3) end position, (4) total length of region:
```
D01-LOCUS100763-290 163 274 111
D01-LOCUS101329-290 37 176 139
D01-LOCUS101812-290 192 298 106
D01-LOCUS104380-290 14 142 128
D01-LOCUS104416-290 154 283 129
D01-LOCUS104617-290 154 298 144
D01-LOCUS10473-290 161 298 137
D01-LOCUS104845-290 27 141 114
D01-LOCUS105006-290 194 298 104
D01-LOCUS10526-290 0 142 142
...
```
