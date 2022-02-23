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

## Stats on locus depth

```
bedtools genomecov -bga -ibam DWR12.final.bam | awk '$4 >= 20' | bedtools merge -i - -c 4 -o min,max,mean,median > DWR12.stats.bed
echo -e "locus\tstart\tend\tmin\tmax\tmean\tmedian\n$(cat DWR12.stats.bed)" > DWR12.stats.bed
```

```
head DWR12.stats.bed 
locus	start	end	min	max	mean	median
D01-LOCUS478-290	14	32	20	21	20.5	20.5
D01-LOCUS478-290	75	87	20	26	23.66666667	24
D01-LOCUS478-290	89	108	24	38	30.85714286	32
D01-LOCUS478-290	109	142	20	40	30.83333333	32
D01-LOCUS478-290	154	258	20	88	70.88461538	78
D01-LOCUS1675-290	0	140	25	188	124.1176471	149
D01-LOCUS2193-290	39	71	27	53	44.16666667	45
D01-LOCUS2193-290	123	143	30	34	32	32
D01-LOCUS2193-290	183	256	21	183	122.7592593	134.5
```


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
