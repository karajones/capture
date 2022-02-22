The scripts and tutorials here are intended to work with the DNA capture bait kit I designed for *Desmognathus* salamanders (family: Plethodontidae). The capture kit itself works with varying degrees of success with other genera I’ve tested it against within Plethodontidae (*Plethodon*, *Eurycea*, and *Gyrinophilus*). The scripts and analyses here should work for any capture kit or similar type of data.  I’ll be adding more scripts and analyses in the near future.

### How the baits were built
The baits were built from a combination of shotgun sequencing and loci from ddRAD. The shotgun sequencing all came from *Desmognathus fuscus*, while the ddRAD loci were primarily *D. quadramaculatus* (northern and southern lineages). The baits are 80 bp long with an overlap of 50% and are intended to map back to reference sequences that are about 300 bp long. The included test data was run on the original kit; I went back afterwards to remove baits that weren’t mapping successfully. The final kit contains:
- 9,756 loci
- 28,256 baits

### Mapping capture reads
[Script](https://github.com/karajones/tutorials/blob/master/scripts/capture_read_mapping.txt) and [tutorial](https://github.com/karajones/tutorials/blob/master/read_mapping.md) for taking raw reads through to clean `bam` files ready for variant calling.
