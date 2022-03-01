# Quality control and read mapping statistics

The ultimate goal here is to curate a list of loci and positions that will produce high quality sets of data for onward analyses. Quality control includes:
1. Making sure SNPs come from a contiguous region of high coverage/depth (not just one small high coverage area or areas with a mix of high/low depth)
2. Removing loci with indels (which are difficult to call accurately) or high numbers of SNPs (see next section)
3. Removing potential duplicate loci/repeating regions

In this section, I'll go over looking at data quality and in the next section about [variant calling](https://github.com/karajones/tutorials/blob/master/vcf_variant_calling.md) I'll cover how to implement quality controls.

## Basic mapping statistics

Outputting the final bam removes all of the information about duplicates, unmapped reads, etc., so if you want to compare the before and after for some basic read mapping stats, then this needs to be run on `.marked.bam` files.
> Note: You can get a lot more detailed output using `stats` rather than `flagstats` but since I'm just mapping back to short sequences rather than trying to align to a chromosome, I find `flagstats` works just fine.
```
for f in *.marked.bam; do samtools flagstat -O tsv $f > ${f%%.*}.flagstats.tsv; done
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
- **primary** - “best” alignment
- **secondary** - alternative alignment (e.g., second “best” alignment if the read maps to more than one place)
- **supplementary** - single read split and aligned to more than one site
- **duplicates** - PCR and optical duplicates
	- “duplicates” includes secondary/supplementary alignments while “primary duplicates” doesn’t
	- I had a lot of duplicates in the example run below because the library was overamplified (which is not ideal!)
- **mapped %** - total reads/mapped
- **primary mapped %** - primary/primary mapped

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

All the secondary, supplementary and duplicates were removed and only primary mapped reads remain, which is good.

# Coverage statistics

>[Bedtools](https://bedtools.readthedocs.io/en/latest/) is required for most of these analyses. See more about Bedtools [genomecov](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html) and [merge](https://bedtools.readthedocs.io/en/latest/content/tools/merge.html). (Those pages have helpful visuals!)

Bedtools `genomecov` outputs statistics on the depth of reads at each individual site on the reference loci. Combined with `merge`, it provides a look at the coverage of reads mapped across each locus. Loci with no mapped reads are removed.
- `genomecov -bg`: output bedgraph format with no zero values included
- `merge -d 500`: regions must be within 500 bp of each other (a number large enough to cover everything on one locus)
- `-c 4 -o min,max,mean,median`: calculate min, max, mean, and median for depth (column 4 on bedgraph output)
- the `echo` command is just used to add a header so the file is easier to read by humans

```
bedtools genomecov -bg -ibam DWR12.final.bam | bedtools merge -d 500 -i - -c 4 -o min,max,mean,median > DWR12.stats.bed &&
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

Example output: (1) genome (rather than individual locus), (2) depth, (3) number of bases at this depth, (4) count of all bases across all loci (will always be the same), (5) fraction of bases at this depth
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

In the above example, there are 1,424,666 bases with zero depth and 247,435 with a depth of two. The highest depth is 28,342. Clearly there are outlier sites with too high of a depth that need to be eliminated. But where to set the limit? Most of the depth falls within a reasonable range. The two dotted lines on the graph represent a 90% and 95% cumulative cut-off; 90% of all bases have a depth less than 14, and 95% of all bases have a depth less than 55.
>I cut the graph off at 200 bp depth so it's easier to see the drop-off of depth.

<img src="https://github.com/karajones/tutorials/blob/master/images/DWR12_depth.png" width="450">

The cut-off points associated with these values can be extracted from the histogram files:

```
awk '{sum += $5; if (sum >= 0.90) {print $2; exit}}' DWR12.coverage.hist.txt
14
awk '{sum += $5; if (sum >= 0.95) {print $2; exit}}' DWR12.coverage.hist.txt
55
```

## Things to consider when looking at depth

Depth of coverage is impacted by three variables that need to be taken into account: sequencing effort, taxonomic relationship, and repeats. Sequencing effort is pretty obvious - the higher the sequencing effort, the higher the total fraction of target bases that have depth and the higher the depth at those bases. Taxa that are more distantly related to the species the capture kit is based on (primarily *Desmognathus fuscus* and *D. quadramaculatus*) will have a higher number of sites with zero coverage, where loci have no aligned reads. This results in a smaller total fraction of target bases covered, regardless of sequencing effort. Here is an example of two taxa at two levels of coverage:

<img src="https://github.com/karajones/tutorials/blob/master/images/depth_by_taxon.png" width="600">

Determining whether a locus is part of a repeating element or transposon based on depth alone is a bit difficult. There are plenty of guides out there that say to just look at, say, average depth and sites with twice (or three times, etc.) are suspect. The logic makes sense intuitively. If baits are capturing from two or more different regions in the genome, those reads will map back to a single locus, producing a locus with twice the depth (or three times or four times, depending on how many different places the sequence shows up. But real data lacks a clear step effect where one locus has regions with depth *x*-times what is found at another locus. This is compounded by the fact that depth tends to be highest in the center of a locus, where more reads overlap, and wane toward the edges.

Long story short, I haven't found a method that is foolproof for dealing with depth. Using a cut-off like the 95% cumulative coverage above makes sense to me as an alternative to just eyeballing a graph of coverage. Other aspects of quality control will catch baits mapping to multiple areas, such as removing reads that map to multiple loci and sites with more than two alleles at a single position.

On to [variant calling](https://github.com/karajones/tutorials/blob/master/vcf_variant_calling.md) >
