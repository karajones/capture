# Mapping capture reads

>This tutorial goes with the [capture_read_mapping.txt](https://github.com/karajones/tutorials/blob/master/scripts/capture_read_mapping.txt) script file.

Each step is set up as a `for` loop to work through all samples before moving on to the next step. Many of the output files are unnecessary - they are set up as checkpoints where extra quality control can be implemented as needed. (Later those files can be sent to stdout and piped on to the next step instead of being saved.) These commands can all be strung together with `&&` to run them consecutively.

### Getting set up
First, install all the programs below. I use [Homebrew] whenever I can if I’m running on a Mac. If you’re on a cluster which doesn’t have these programs set up as modules, links to the binaries are below:
- `Trimmomatic`: [https://github.com/usadellab/Trimmomatic]
- `BWA`: [https://github.com/lh3/bwa]
- `samtools`: [https://github.com/samtools/samtools]
- `picard-tools`: [https://broadinstitute.github.io/picard/]
>Note: I’m not using BWA-MEM2 here because I’m using an M1 Mac and it’s not working with the new chip architecture yet. BWA installs with an M1-specific following changes described here: [https://stackoverflow.com/questions/65555170/unable-to-run-make-command-for-bwa-on-apple-m1-on-mac-os-big-sur]

### Initial sequence quality control
Trim and remove adapters using Trimmomatic:
- `ILLUMINACLIP:TruSeq3-PE.fa:2:30:10`: remove adapters
- `2:True`: keep both paired reads
- `MINLEN:40`: drop reads shorter than 40 bases long
```
for f in ./fastqs/*_1.fq.gz; do trimmomatic PE ${f%%_*}_1.fq.gz ${f%%_*}_2.fq.gz ${f%%_*}_1p.fq.gz ${f%%_*}_1u.fq.gz ${f%%_*}_2p.fq.gz ${f%%_*}_2u.fq.gz ILLUMINACLIP:./adapters/TruSeq3-PE.fa:2:30:10:2:True MINLEN:40; done
```

### 0. Set up new directories
All subdirectories are assumed to reside within one main directory. To run the code as-is, within one directory you should have:
- Raw files are stored in `fastqs`
- Indexed reference fasta stored in `reference`
- Trimmomatic adapters stored in `adapters`
```
mkdir sams
mkdir bams
```

### 1. Align sequences to reference using BWA MEM
If you want to use a different aligner, keep in mind that the quality control values used below will change since they are specific to bwa mem.
> Note: I really like files neatly stored and named, so I make heavy use of `#` and `%` operators in bash/zsh to strip extensions and make prettier filenames. 
```
for f in ./fastqs/*_1p.fq.gz; do bwa mem ./reference/reference.fas $f ${f%%_*}_2p.fq.gz > ./sams/${${f%%_*}##*/}.sam; done
```

### 2. Convert to bam (-b) and sort
I don’t remove any unmapped sequences at this point so I can run statistics on the bams below and actually get info on mapped vs unmapped. If you want to remove unmapped reads to reduce file size and time running MarkDuplicates add `-F 4` to `samtools view` call. 

> NOTE: mapq is confusing for bwa-mem (there are no good values that will just remove "low quality" mappings/reads) so I don't specify a -q value here.
- `samtools view`: 
	- `-h` include header (needed for next steps)
	- `-u` uncompressed bam (so samtools doesn’t have to compress and uncompress before the next step)
	- `-b` output bam
- `samtools sort`: `-m` increase memory to 8G
```
for f in ./sams/*.sam; do samtools view -h -u -b $f | samtools sort -m 8G - > ./bams/${${f%.*}##*/}.bam; done
```

### 3. Mark duplicates with Picard
The total length of an aligned read can be longer than the template sequence because of overhang beyond the original sequence length. `--MINIMUM_DISTANCE 500`  tells Picard to search for duplicates within a window of 500 bp rather than the default of `-1` which is limited to the first read’s read length. In other words, if `--MINIMUM_DISTANCE -1` and the first read is 150 bp long, then the window to search for duplicates within an alignment is 300 bp, which isn’t long enough to encompass the maximum alignment length.
> Increasing the minimum distance also increases the amount of memory used but shouldn’t cause problems on any decent computer.
```
for f in ./bams/*.bam; do picard MarkDuplicatesWithMateCigar --MINIMUM_DISTANCE 500 -I $f -O ./bams/${${f%.*}##*/}.marked.bam -M  marked_duplicates.txt; done
```

### 4. Output clean bam file
Remove all the unwanted reads from the alignment and output a `.final.bam` file. Some of the calls used here are probably redundant but better safe than sorry.
- `samtools view`: `-h` include header 
- `-F 1284` flags used () to exclude reads: 
	- `0x4` read unmapped
	- `0x100` not primary alignment
	- `0x400` read is PCR or optical duplicate. 
>See helpful page on samtools flags here: [https://broadinstitute.github.io/picard/explain-flags.html]

The `grep -e` call excludes alternative alignments (XA) and split alignments (SA). I don’t want reads that map to more than one location because I have no way to check whether they’re valid or not since there’s no reference genome. (The reference are just random ddRAD and shotgun sequences so they could be parts of repeating elements or other ephemera.) Flags only remove the secondary alignment, not the original alignment so `grep` is needed to remove the read completely from all locations where it was mapped.  Apparently this can also be done using `sambamba` but I didn’t bother to download another program and this works just fine. `¯\_(ツ)_/¯`
> For more Info, see this discussion: [https://bioinformatics.stackexchange.com/questions/508/obtaining-uniquely-mapped-reads-from-bwa-mem-alignment]
```
for f in ./bams/*.marked.bam; do samtools view -h -F 1284 $f | grep -v -e 'XA:Z:' -e 'SA:Z:' - | samtools view -b - > ./bams/${${${f%.*}%.*}##*/}.final.bam; done && for f in ./bams/*.final.bam; do samtools index $f; done
```

### 5.  Index sorted bam files
This is needed for some onward steps.
```
for f in *.final.bam; do samtools index $f; done
```
