# Variant calling and VCF output

### Create an individual VCF for each sample
1. `mpileup -d`: maximum depth; see the [quality control and statistics](https://github.com/karajones/tutorials/blob/master/quality_control_statistics.md) section
2. `call`: call variants
   - `-m` use multi-allelic caller (recommended)
   - `-v` output only variant sites
3. `view`
   - `-m2` and `-M2` set min and max number of alleles to 2
4. `filter`
   - `--IndelGap 5`: assuming indels within 5 bp of each other are false. 
     - See: [https:/ ls.github.io/bcftools/howtos/consensus-sequence.html](https://samtools.github.io/bcftools/howtos/consensus-sequence.html)
   - `-e`: exclude calls with a quality score below 200 and a depth below 20
5. `norm`: normalize indels
   - see explanation here: [https://genome.sph.umich.edu/wiki/Variant_Normalization](https://genome.sph.umich.edu/wiki/Variant_Normalization)
   - `-d all` in the case of records with the same position, only the first will be considered and appear on output.
   - `-c wx` warn if REF allele is incorrect or missing and exclude

```
for f in ./bams/*.final.bam; do bcftools mpileup -d 500 -f ./reference/reference_updated.fas $f | bcftools call -m -v | bcftools view -m2 -M2 | bcftools filter -e '%QUAL<200 || DP<20' | bcftools norm -f ./reference/reference.fas - -c wx > ${f%%.*}.20dp.vcf; done
```

Example log (normal):
```
Note: none of --samples-file, --ploidy or --ploidy-file given, assuming all sites are diploid
[mpileup] 1 samples in 1 input files
[mpileup] maximum number of reads per input file set to -d 874
NON_ACGTN_REF	D07-LOCUS145533-290	241	WG
NON_ACGTN_REF	D01-LOCUS40720-290	62	GCRGGAAA
NON_ACGTN_REF	GA03F-LOCUS15439-135	91	YGGTCTGGAAT
NON_ACGTN_REF	GA06F-LOCUS91669-139	80	TAYTACA
NON_ACGTN_REF	GA03F-LOCUS134434-135	60	GGATYTTAGGGTAGGTA
NON_ACGTN_REF	DF-JK13-1-010-LOCUS60855-290	130	AGSCCCCCTCCNNNNNNNNNNNNNNCMG
NON_ACGTN_REF	DC-JK13-112-LOCUS61301-290	125	CRGTTTTATTACAATCANNNNNNNNNNNNRGTTTTATTACAATCA
NON_ACGTN_REF	DF-JK13-1-010-LOCUS97135-290	217	GTCCTRGAGA
NON_ACGTN_REF	DF-JK13-1-010-LOCUS113290-290	138	YAACNNNNNNNNNNNNTTTAAACTAYAAC
NON_ACGTN_REF	DF-JK13-1-012-LOCUS116598-290	134	ACRTTGGTNNNNNNNNNNNNCRTTGGT
NON_ACGTN_REF	DC-JK13-112-LOCUS146502-290	138	CRGGNNN
NON_ACGTN_REF	DC-JK13-112-LOCUS156198-290	129	GGCYGGATAAATACNNNNNNNNNNNNNGCYGGATAAATAC
NON_ACGTN_REF	DC-JK13-130-LOCUS18579-290	53	TTCRT
Lines   total/split/realigned/skipped:	125975/0/1292/13
```

### Create a multi-sample VCF
Calls as explained above but creating a VCF from all `bam` files in a directory.

```
for f in *.final.bam; do echo $f; done > bamlist &&
bcftools mpileup -d 500 -f ./reference/reference_updated.fas $f | bcftools call -m -v | bcftools view -m2 -M2 | bcftools filter -e '%QUAL<200 || DP<20' > all_samples.vcf
```

## Subsetting the VCF using whitelists and blacklists 
The VCF output above has basic quality controls, but extra quality measures can be implemented by searching for loci that match a set of conditions, creating a blacklist of sites to exclude (or a whitelist of sites to include), and then using that list to subset the original VCF. `VCFtools` is really handy for this and can include/exclude sites for an entire locus (`--chr`/`--not-chr`), a range of positions on a locus (`--from-bp` and `--to-bp`), and exact position (`--positions`/`exclude-positions`). The ultimate goal would be to create a whitelist of loci and corresponding positions to pull out the same SNP sites from each individual every time. 

### Find contiguous regions of high coverage/depth
Here's an example of how to create a whitelist of sites. I use `bedtools genomecov` on the final bam file to find all regions with non-zero coverage, `awk` to select only regions with at least 50 bp depth, and the `merge` function to combine any of these high depth regions of coverage that abut or overlap. Then I use awk to check for several conditions:
1. `$3-$2>=100`: length of the region must be greater than or equal to 100 contiguous bp
2. `$6<400`: mean read mapping depth must be less than 400 bp
3. `$2==0`: regions must start at 0 (mapped to the beginning of a locus
> `awk` is acting on the columns from the `merge` output: (1) locus id, (2) start position, (3) end position, (4) minimum depth, (5) maximum depth, (6) mean depth. See the stats section.

```
bedtools genomecov -bga -ibam H1.final.bam | sort -k1,1 -k2,2n | awk '$4 >= 50' | bedtools merge -i - -c 4 -o min,max,mean | awk -v OFS='\t' '{if($3-$2>=100 && $6<400 && $2==0) print $1,$2,$3}' > whitelist.bed
```
Output columns: (1) locus name, (2) start position, (3) end position
```
D01-LOCUS11635-290	0	142
D01-LOCUS116965-290	0	105
D01-LOCUS12397-290	0	125
D01-LOCUS127494-290	0	144
D01-LOCUS131414-290	0	141
D01-LOCUS131614-290	0	132
D01-LOCUS136623-290	0	142
D01-LOCUS139306-290	0	132
D01-LOCUS147682-290	0	142
D01-LOCUS166291-290	0	119
...
```

Use the whitelist in `VCFtools` to create file with subset of SNPs:
>`VCFtools` will always output the new file as `name.recode.vcf`. So below, the output file would be `whitelisted_variants.recode.vcf`.
```
vcftools --vcf desmog_variants.vcf --positions whitelist.bed --recode --out whitelisted_variants
```

### Limit to one region per locus
Uses the `bed` file produced above to create a file which lists *one* region for each locus where all bp in a region are at least 100 bp long and 50 bp depth (or whatever was specified in the `bed` file). The file is sorted (`sort`) on locus id (`-k1,1`) and then only unique (`-u`) locus id's are retained. This command can be piped with the previous command to output this file without the intermediate `bed` file. The resulting file can be used the same as the whitelist as shown above.

```
cat whitelist.bed | sort -u -k1,1 > whitelist_unique.bed
```

## Subset one random SNP per locus from VCF
Some analyses (like those using a site frequency spectrum) require all SNPs, but most SNP-based analyses require unlinked SNPs. This command will print a whitelist of loci/position with one *random* SNP per locus:

>NOTE: Macs and some flavors of Linux might not have the `shuf` command. If not, install `coreutils`. (i.e., `brew install coreutils`) Every time the command is run the output will be a different set of locus/position combinations!

1. pull out all loci and positions from VCF
2. `shuf` randomly shuffle the output from the first awk command
3. outputs first unique value from column 1 (locus ID)
```
sed -e 's/chr//' H1.vcf | awk -F"\t" '{if (!/^#/) {print $1,$2}}' | shuf | awk '!a[$1]++' > H1.random.snps
```

Output (locus, position):

```
GA04F-LOCUS173676-135 44
D07-LOCUS100541-290 35
GA03F-LOCUS112744-135 27
GA03F-LOCUS36376-135 120
GA05F-LOCUS38950-135 50
GA03F-LOCUS160728-135 68
D01-LOCUS94303-290 44
GA03F-LOCUS46506-135 114
D07-LOCUS100261-290 135
GA04F-LOCUS95349-135 36
...
```

## Print list of loci that have exactly 3 SNVs (or something similar)
This can easily be modified to output different levels of SNVs. For example, `if(a[i]<10)` will output all loci with fewer than 10 SNVs.

```
sed -e 's/chr//' all_variants.vcf | awk -F"\t" '{if (!/^#/) a[$1]++}END{for(i in a) if(a[i]=3) print i}' > 3snps.txt
```

Output (locus):

```
GA03F-LOCUS44720-135
DC-JK13-112-LOCUS69498-290
DF-JK13-1-010-LOCUS95860-290
DF-JK13-1-015-LOCUS93698-290
D13-LOCUS107014-290
DF-JK13-1-010-LOCUS122818-290
D07-LOCUS55782-290
GA03F-LOCUS160728-135
D07-LOCUS35049-290
D07-LOCUS120071-290
...
```

Since the output is a list of loci, these can be whitelisted in `VCFtools` using the `--chr` flag (or the `--not-chr` flag if you want to exclude the loci):
```
vcftools --vcf desmog_variants.vcf --chr 3snps.txt --recode --out 3snps
```

## List loci with indels
Print list of loci that contain indels (output can be used as a blacklist to remove indels later).
>Indels can be removed in VCFtools with the `--remove-indels` flag (and with `bcftools` during variant calling), but that only removes the indels themselves, not the entire locus that an indel is found on. Since indels are difficult to call and throw off the positioning of other variant sites on a locus, it may be wise to completely remove loci with indels.

```
sed -e 's/chr//' all_samples.vcf | awk '{OFS="\t"; if (!/^#/ && /INDEL/){print $1}}' | uniq > blacklist
```

Output (locus):
```
D08-LOCUS50108-290
D07-LOCUS54179-290
D07-LOCUS109722-290
D07-LOCUS120960-290
D07-LOCUS129094-290
D08-LOCUS133320-290
D08-LOCUS133882-290
D01-LOCUS143252-290
D07-LOCUS148584-290
D07-LOCUS162136-290
```

This will output a VCF with all loci that contain indels removed:
```
vcftools --vcf desmog_variants.vcf --chr blacklist --recode --out no_indels
```

## Various VCF statistics

These commands can be used for either individual or multi-sample VCFs.

### Print list of number of SNVs (SNPs and indels) per locus

```
sed -e 's/chr//' all_variants.vcf | awk '{if (!/^#/) a[$1]++}END{for(i in a) print i,a[i]}'
```

Example output (locus, number of SNVs):
```
GA03F-LOCUS44720-135 3
DC-JK13-112-LOCUS69498-290 15
DF-JK13-1-010-LOCUS95860-290 4
DF-JK13-1-015-LOCUS93698-290 14
D13-LOCUS107014-290 4
DF-JK13-1-010-LOCUS122818-290 2
D07-LOCUS55782-290 5
GA03F-LOCUS160728-135 1
D07-LOCUS35049-290 4
D07-LOCUS120071-290 17
...
```

### Print SNV count at each locus

```
sed -e 's/chr//' wrighti_all_variants.vcf | awk '{if (!/^#/) a[$1]++}END{for(i in a) print a[i]}' | sort -n | uniq -c
```

Output (number of loci, number of SNPs):
```
173 1
  94 2
  89 3
  83 4
  59 5
  70 6
  48 7
  32 8
  32 9
  22 10
  42 11
  16 12
...
```

### Output proportion of missing data
For multi-sample VCFs only. Output the proportion of missing data for each individual/sample listed in the VCF.

```
paste \
<(bcftools query -f '[%SAMPLE\t]\n' all_variants.vcf | head -1 | tr '\t' '\n') \
<(bcftools query -f '[%GT\t]\n' all_variants.recode.vcf | awk -v OFS="\t" '{for (i=1;i<=NF;i++) if ($i == "./.") sum[i]+=1 } END {for (i in sum) print i, sum[i] / NR }' | sort -k1,1n | cut -f 2)
```




