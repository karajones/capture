> [!WARNING]
> Manuscript in review. Some of this might be out of date - I'll update this when I have time.

# Tutorials and scripts for plethodontid capture kit

The scripts and tutorials here are intended to work with the DNA capture bait kit I designed for *Desmognathus* salamanders (family: Plethodontidae). The capture kit itself works with varying degrees of success with other genera I’ve tested it against within Plethodontidae (*Plethodon*, *Eurycea*, and *Gyrinophilus*). The scripts and analyses here should work for any capture kit or similar type of data.  I’ll be adding more scripts and analyses in the near future.

### How the baits were built
The baits were built from using loci from ddRAD, primarily from *Desmognathus fuscus* and *D. quadramaculatus* (northern and southern lineages). The baits are 80 bp long with an overlap of 50% and are intended to map back to reference sequences that are about 300 bp long. The included test data was run on the original kit; I went back afterwards to remove baits that weren’t mapping successfully. The final kit contains:
- 9,756 loci
- 28,256 baits

### Lab work
[Sort-of tutorial](https://github.com/karajones/capture/blob/0bc35a31d067d0059ff2feafffae87a8fbd8b947/lab_work.md) on how to get raw DNA ready for sequencing. Includes list of reagents and recommendations for changes to standard protocols.

### Mapping capture reads
[Script](https://github.com/karajones/tutorials/blob/master/scripts/capture_read_mapping.txt) and [tutorial](https://github.com/karajones/tutorials/blob/master/read_mapping.md) for taking raw reads through to clean `bam` files ready for variant calling.

### Quality control and read mapping statistics
[Tutorial](https://github.com/karajones/tutorials/blob/master/quality_control_statistics.md) for producing read mapping statistics and thinking about quality control measures to implement.

### Variant calling and VCF output
[Script](https://github.com/karajones/tutorials/blob/master/scripts/vcf_script.txt) and [tutorial](https://github.com/karajones/tutorials/blob/master/vcf_variant_calling.md) for producing a single or multi-sample VCF. How to use whitelist, blacklists, and other tools to implement further quality control measures on the output. Also, how to get quick statistics from a VCF.

### Automated fastq to VCF shell script
[Script](https://github.com/karajones/capture/blob/master/scripts/capture_multi_vcf.sh) starts with a directory of trimmed (paired) fastqs and outputs multi-sample VCF for all the samples in the directory with one random SNP per locus. Determines max depth based on average depth across all sites, removes loci with indels, keeps only sites where 50% or more of individuals have data. Requires: bwa, samtools, picard, bedtools, bcftools, vcftools.

### To do:
- Intersecting bed files to create a whitelist of loci and positions
- Output a fasta file with variants for phylogenetic analyses
- Check for linkage between loci/finish mapping to transcriptome
- Create single output file for mapping statistics (from individual flagstat files)
- Optimize for parallelizing on computing cluster
