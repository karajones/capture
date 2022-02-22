# Variant calling: VCF

### Create an individual VCF for each sample
1. `mpileup -d`: maximum depth; see stats section
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
for f in ./bams/*.final.bam; do bcftools mpileup -d 874 -f ./reference/reference_updated.fas $f | bcftools call -m -v | bcftools view -m2 -M2 | bcftools filter -e '%QUAL<200 || DP<20' | bcftools norm -f ./reference/reference.fas - -c wx > ${f%%.*}.20dp.vcf; done
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
bcftools mpileup -d 874 -f ./reference/reference_updated.fas $f | bcftools call -m -v | bcftools view -m2 -M2 | bcftools filter -e '%QUAL<200 || DP<20' > all_samples.vcf
```

## VCF statistics

These commands can be used for either the individual or multi-sample VCF.

### Print list of number of SNVs (SNPs and indels) per locus:

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

### Print SNV counts:

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


### Indels
Print list of loci that contain indels (these will be blacklisted):

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

## Subset one random SNP per locus from VCF
Print a whitelist of loci/position with one *random* SNP per locus:

>NOTE: Macs and some flavors of Linux might not have the `shuf` command. If not, install `coreutils`. (i.e., `brew install coreutils`) Every time the command is run the output will be a different set of locus/position combinations!

1. pull out all loci and positions from VCF
2. `shuf` randomly shuffle the output from the first awk command
3. outputs first unique value from column 1 (locus ID)
```
sed -e 's/chr//' H1.recode.vcf | awk -F"\t" '{if (!/^#/) {print $1,$2}}' | shuf | awk '!a[$1]++' > H1.random.snps
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

Print list of loci that have exactly 3 SNVs:
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
