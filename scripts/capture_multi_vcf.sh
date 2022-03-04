#!/bin/sh

# NOTE: required programs/modules
# bwa, samtools, picard, bedtools, bcftools, vcftools

# NOTE: the reference should already be indexed
reference=./reference/reference_updated.fas

# make all the directories
mkdir -p sams bams tmp final_bams stats

# create file for depths
depths=./tmp/depths
touch $depths

# create bamlist
bamlist=./tmp/bamlist
touch $bamlist

# NOTE: trimmed fastqs assumed to be in a directory name `trimmed` (with .fq.gz extensions)
for f in ./trimmed/*_R1_paired.fq.gz
    do

    base=$(basename ${f%%_*})
    R1=./trimmed/${base}_R1_paired.fq.gz
    R2=./trimmed/${base}_R2_paired.fq.gz
    sam=./sams/${base}.sam
    bam=./bams/${base}.bam
    marked=./bams/${base}.marked.bam
    dups=./tmp/${base}_marked_duplicates.tmp
    final=./final_bams/${base}.final.bam
    flags=./stats/${base}.flagstats.tsv
    whitelist=./tmp/${base}.whitelist
    blacklist=./tmp/${base}.blacklist
    random=./tmp/${base}.random.snps
    ID="ID:"${base}
    SM="SM:"${base}

    # add file to list of samples for VCF
    echo $final >> $bamlist

    echo "(>‘o’)> Mapping $base to reference..."
    bwa mem $reference $R1 $R2 > $sam
    
    echo "<('.')> Coverting to bam and sorting..."
    samtools view -h -u -b $sam | samtools sort -m 8G - > $bam

    echo "^('-')^ Marking duplicates with picard..."
    picard MarkDuplicatesWithMateCigar --MINIMUM_DISTANCE 500 -I $bam -O $marked -M $dups

    echo "<(｀^´)> Outputting final bam..."
    # addreplacerg is needed for nicer multi-sample vcf sample names
    samtools view -h -F 1284 $marked | grep -v -e 'XA:Z:' -e 'SA:Z:' - | samtools addreplacerg -r $ID -r $SM - | samtools view -b - > $final

    # index final bam
    samtools index $final

    echo "(^-^*)/ Outputting flag statistics file..."
    samtools flagstat -O tsv $marked > $flags

    # calculate high depth cut-off
    depth=$(bedtools genomecov -ibam $final | grep ^genome | awk '{sum += $5; if (sum >= 0.99) {print $2; exit}}')
    echo "(-_q) 99% of sites have a depth less than $depth"

    # write depth to file for later use
    echo $depth >> $depths 

    echo "Processing for $base is finished! =^_^= \n"

done

# output multi-sample vcf

# vcf variables
vcf=all_samples.vcf
final_vcf=all_samples

# calculate average depth across individuals
depth2=$(cat $depths | awk '{sum += $1} END {print sum/NR}')
echo "( •_•) Average depth across individuals is $depth2"

# output multi-sample vcf
echo "(*'ｰ')/ Creating multi-sample VCF..."
bcftools mpileup -d $depth2 -f $reference -b $bamlist | bcftools call -m -v | bcftools view -m2 -M2 | bcftools filter -e '%QUAL<200 || DP<20' > $vcf

# create list of loci with indels
echo ">^_^< Applying whitelist to VCF..."
sed -e 's/chr//' $vcf | awk '{if (!/^#/ && /INDEL/){print $1}}' | uniq > $blacklist

# create list of random snps
sed -e 's/chr//' $vcf | awk -F"\t" '{if (!/^#/) {print $1,$2}}' | shuf | awk '!a[$1]++' > $random 

#remove loci with indels from whitelist
awk 'NR==FNR{a[$0];next} !($1 in a)' $blacklist $random > $whitelist

# output final vcf
# max-missing: keep only sites where 50% or more of individuals have data
vcftools --vcf $vcf --positions $whitelist --max-missing 0.5 --recode --out $final_vcf

# remove temporary files
rm -rf ./tmp 
rm -f $vcf

echo "＼(^o^)／ Finished! \nThe final VCF is called $final_vcf.recode.vcf"
