# tssv-library
A small tool to generate tssv library files from a VCF file.

## Installation
This assumes you have Conda installed.

```bash
# Create the conda environment
conda create -n tssv-library
# Activate the conda environment
conda activate tssv-library

# Install the dependencies
conda install pysam cyvcf2
```

## Usage
tssv-library requires a VCF file, and an indexed reference file as input.

### Index the reference file
tssv-library requires a fasta file indexed with samtools as input.
```bash
samtool index reference.fasta
```

### Options
```
usage: create-library.py [-h] --reference REFERENCE --vcf VCF [--flank-size FLANK_SIZE] [--max-size MAX_SIZE]

optional arguments:
  -h, --help            show this help message and exit
  --reference REFERENCE
                        Reference fasta file (default: None)
  --vcf VCF             Create a library for the variants in the VCF (default: None)
  --flank-size FLANK_SIZE
                        Size of the flanking regions (default: 20)
  --max-size MAX_SIZE   Maximum size of an indel to verify (default: 20)
```

### Caveats
Make sure that you do not request a library that is too large for your read
size. If `2*(--flank-size) + --max-size` is larger than the read size of your
sequencing library, you will never find any hits. Also note that if this is
only slightly smaller than your read size, very few reads will perfectly
overlap the requested region, so you might still get very few results.
