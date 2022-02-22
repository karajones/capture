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
