#!/usr/bin/python
# coding=UTF-8

import click
import re
import os
import math
import numpy as np
import gzip

@click.command()
@click.option("--vfile","-v",help="variation file in special format")
@click.option("--rfile","-r",help="reference file in fasta format, hg38")

def main(rfile, vfile):
    print('Gene\tChrom\tPos\tRef\tAlt\tMutSeq\tWholeSeq\tRegion\tStrand\tTranscript\tUTR-len')
    GeneInf, NM = {}, {}
    pwdir = os.path.dirname(os.path.abspath(__file__))
    mfile = pwdir+'/database/MANE/MANE.GRCh38.v1.1.ensembl_genomic.gff.gz'
    with gzip.open(mfile) as mf:
        for mline in mf:
            mline = str(mline, encoding='utf-8')
            if mline.startswith('#'):
                continue
            minf = mline.strip().split('\t')
            mchr, mfunc, mstart, mend, mstrand, mdetail = minf[0], minf[2], minf[3], minf[4], minf[6], minf[8]
            if 'UTR' not in mfunc and 'transcript' not in mfunc:
                continue
            mgene, menst = mdetail.split('gene_name=')[1].split(';')[0], mdetail.split('transcript_id=')[1].split('.')[0]
            if 'transcript' in mfunc:
                mnmid = mdetail.split('RefSeq:')[1].split('.')[0]
                NM[menst] = mnmid
            GeneInf.setdefault(mgene, {})
            GeneInf[mgene]['base'] = mchr+'\t'+mstrand
            GeneInf[mgene].setdefault(menst, {})
            GeneInf[mgene][menst].setdefault('UTR5', [])
            GeneInf[mgene][menst].setdefault('UTR3', [])
            if mfunc == 'five_prime_UTR':
                GeneInf[mgene][menst]['UTR5'].append(mstart+'\t'+mend)
            elif mfunc == 'three_prime_UTR':
                GeneInf[mgene][menst]['UTR3'].append(mstart+'\t'+mend)
                            
    reffa = rfile
    with open(vfile) as vf:
        for vline in vf:
            if vline.startswith('var'):
                continue
            vinf = vline.strip().split('\t')
            vchr, vpos, vref, valt, vgene, vregion = vinf[0], vinf[1], vinf[2], vinf[3], vinf[4], vinf[5]
            if vgene not in GeneInf:
                continue
            chrom, strand = GeneInf[vgene]['base'].split('\t')[0], GeneInf[vgene]['base'].split('\t')[1]
            utr, utr_seq = 'UTR5', {}
            if '3_prime_UTR_variant' in vregion:
                utr = 'UTR3'
            for enst in GeneInf[vgene]:
                if enst == 'base':
                    continue
                transcript = enst+':'+NM[enst]
                utr_inf = GeneInf[vgene][enst][utr]
                for i in range(0, len(utr_inf)):
                    utr_s, utr_e = utr_inf[i].split('\t')[0], utr_inf[i].split('\t')[1]
                    if int(vpos) < int(utr_s) or int(vpos) > int(utr_e):
                        continue
                    region = chrom+':'+utr_s+'-'+utr_e
                    seq = os.popen("samtools faidx %s %s" % (reffa, region)).read().strip().split('\n')
                    seq = ''.join(seq[1:len(seq)])
                    if strand == '-':
                        seq = revcomp(seq)
                    ulen = str(int(utr_e) - int(utr_s) + 1)
                    mutseq = Mutseq_get(utr_s, utr_e, seq, vpos, vref, valt, strand)
                    #print(utr_s+'\t'+utr_e)
                    print(vgene+'\t'+chrom+'\t'+vpos+'\t'+vref+'\t'+valt+'\t'+mutseq+'\t'+seq+'\t'+vregion+'\t'+strand+'\t'+transcript+'\t'+ulen)

def Mutseq_get(utr_s, utr_e, seq, vpos, vref, valt, strand):
    p = int(vpos)-int(utr_s)
    u_seq, d_seq = seq[0:p], seq[p+len(vref):len(seq)]
    if strand == '-':
        p = int(utr_e)-int(vpos)
        u_seq, d_seq = seq[0:p-len(vref)+1], seq[p+1:len(seq)]
        vref, valt = revcomp(vref), revcomp(valt)
    mseq = u_seq+'('+vref+'/'+valt+')'+d_seq
    return mseq

def complement(seq2):
    trantab = seq2.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh')
    return seq2.translate(trantab)

def revcomp(seq1):
    return complement(seq1)[::-1]
                    

if __name__ == "__main__":
    main()
