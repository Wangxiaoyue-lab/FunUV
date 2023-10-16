# FunUV
A prediction tool for functional UTR variants including UTR variation annotation.

# Dependencies
Python 3.9 and standard packages (scipy, sklearn, numpy, pandas, pyBigWig). 
SAMTools 1.16
ViennaRNA-2.4.9
MEME-5.3.3

# Installation
Clone this github repository.

# Usage
1. annotate variation using region-specific features.

  input file format (tab split):

  Chrom   Pos     Ref     Alt     Gene    Region
  chr1    1013541 T       C       ISG15   5_prime_UTR_variant
  chr1    1232241 G       A       B3GALT6 5_prime_UTR_variant

python3 Sub_MutSeqByMANE.py -v inputfile -r Homo_sapiens_assembly38.fasta

2.annotate variation using region-specific features.
