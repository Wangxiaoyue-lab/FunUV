# FunUV
A prediction tool for functional UTR variants including UTR variation annotation.

# Dependencies
Python 3.9 and standard packages (scipy==1.10.1, scikit-learn==1.2.1, numpy==1.23.5, pandas==1.1.4, pyBigWig==0.3.21).  
SAMTools-1.16  
ViennaRNA-2.4.9  
MEME-5.3.3  

# Installation
Clone this github repository.  
The relevant conservation score files are too large to upload, please downloaded from [hg38.phastCons100way.bw](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.phastCons100way.bw) and [hg38.phyloP100way.bw](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw).  Then place the two files in the Database/Conservation directory.  
The prediction model for UTR3 is too large, please get it from the release and then place in the Model directory  

Note: All analysis is based on the hg38 version of the reference genome.  


# Usage
1. prepare transcript information.  

    input file format (tab split):  

    Chrom  Pos  Ref  Alt  Gene  Region  
    chr1  1013541  T  C  ISG15  5_prime_UTR_variant  
    chr1  1232241  G  A  B3GALT6  5_prime_UTR_variant

    python3 Sub_MutSeqByMANE.py -v inputfile -r Homo_sapiens_assembly38.fasta > test_5utr.addseq.txt

    output file format (tab split):

    Gene    Chrom    Pos    Ref    Alt    MutSeq    WholeSeq    Region    Strand    Transcript    UTR-len
    ISG15    chr1    1013541    T    C    GGCGGCTGAGAGGCAGCGAACTCATCTTTGCCAGTACAGGAGCT(T/C)GTGCCGTGGCCCACAGCCCACAGCCCACAGCC    GGCGGCTGAGAGGCAGCGAACTCATCTTTGCCAGTACAGGAGCTTGTGCCGTGGCCCACAGCCCACAGCCCACAGCC    5_prime_UTR_variant    +    ENST00000649529:NM_005101    77  
    B3GALT6    chr1    1232241    G    A    ACTC(G/A)CGAGTCCGGCCTGGGCCGCCGGCCCGGCGCGGGCGCC    ACTCGCGAGTCCGGCCTGGGCCGCCGGCCCGGCGCGGGCGCC    5_prime_UTR_variant    +    ENST00000379198:NM_080605    42  

2. annotate variant using region-specific features.  

    input file: test_5utr.addseq.txt  
    
    python3 FunUV_annotation.py -f test_5utr.addseq.txt -o test_5utr -r 5UTR  
    
    output file format (5UTR):  

    Var     Seq     PCS     PPS     MFE     IRES    uAUG    uStop   uORF    Kozak-score     GC-ratiAC-ratio AG-ratio        AT-ratio        TC-ratio        TG-ratio        A-ratio T-ratio C-ratioG-ratio  AA-count        AT-count        AC-count        AG-count        TA-count        TT-counTC-count TG-count        CA-count        CT-count        CC-count        CG-count        GA-counGT-count GC-count        GG-count        A-homo  T-homo  C-homo  G-homo  AA-dimo AT-dimo AC-dimoAG-dimo  TA-dimo TT-dimo TC-dimo TG-dimo CA-dimo CT-dimo CC-dimo CG-dimo GA-dimo GT-dimo GC-dimoGG-dimo  sequni  
    ISG15_chr1_1013541_T_C_5UTR_44  ACTCATCTTTGCCAGTACAGGAGCT(T/C)GTGCCGTGGCCCACAGCCCACAGC  0.0    -0.12399999797344208     1.5999999999999996      0       0       0       0.0     0       0.62   0.58     0.44    0.38    0.56    0.42    0.2     0.18    0.38    0.24    0       1       4      10  
    B3GALT6_chr1_1232241_G_A_5UTR_4 ACTC(G/A)CGAGTCCGGCCTGGGCCGCCGGCC       0.0     0.026000000536441803    1.0999999999999996      0       0       0       0.0     0       0.7931034482758621     0.5517241379310345       0.4482758620689655      0.20689655172413793     0.5517241379310345     0.4482758620689655       0.10344827586206896     0.10344827586206896     0.4482758620689655     0.3448275862068966       0       0       2       1       0       0       2       1       1  

    output file format (3UTR):  

    Var     Seq     PCS     PPS     MFE     miRNA   PAS     AU-cls  AU-9mer AU-13mer        GU     CU       Pumilio AU      AU-elements     GC-ratio        AC-ratio        AG-ratio        AT-ratiTC-ratio TG-ratio        A-ratio T-ratio C-ratio G-ratio AA-count        AT-count        AC-counAG-count TA-count        TT-count        TC-count        TG-count        CA-count        CT-counCC-count CG-count        GA-count        GT-count        GC-count        GG-count        A-homo T-homo   C-homo  G-homo  AA-dimo AT-dimo AC-dimo AG-dimo TA-dimo TT-dimo TC-dimo TG-dimo CA-dimoCT-dimo  CC-dimo CG-dimo GA-dimo GT-dimo GC-dimo GG-dimo sequni
    AGRN_chr1_1055037_T_C_3UTR_55   CACCAGAGCCCCGCGCCCGCTGTAATTATTTTCTATTTTTGTAAACTTGT(T/C)GCTTTTTGATATGATTTTCTTGCCTGAGTGTTGGCCGGAGGGACTGCTG        1.0     2.302999973297119       2.3999999999999986      0       0       0       0       0       0       0       0       0       0       0.47   0.39     0.4     0.53    0.6     0.61    0.16    0.37    0.23    0.24    3       6       3      18       3       11      2       8       8       5       6       5       8       4       3      33
    AGRN_chr1_1055153_G_A_3UTR_171  GTCCAGGCAGCCGTGCTGCAGACAGACCTAGTGCCGAGGGATGGACAGGC(G/A)AGGTGGCAGCGTGGAGGGCTCGGCGTGGATGGCAGCCTCAGGACACACA        0.0     -2.4809999465942383     0.7999999999999972      0       0       0       0       0       0       0       0       0       0       0.66   0.48     0.62    0.34    0.38    0.52    0.22    0.12    0.26    0.4     1       2       6      12       1       0       3       8       12      4       5       5       8       6       12     14       2       1       2       3       1       1       3       1       1       0       1      20  
3. predict variant effect using region-specific features.

   input file: test_5utr_annotation.txt  

   python3 Prepare_forPred.py -v test_5utr_annotation.txt -r 5UTR > test_5utr_subanno.txt
   python3 FunUV_predict.py -pf test_5utr_subanno.txt -r 5UTR -o test_5utr
