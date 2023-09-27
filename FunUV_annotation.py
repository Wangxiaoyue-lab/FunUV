#!/usr/bin/python

import click
import re
import os
import math
import numpy as np
import pyBigWig
from sklearn.preprocessing import OneHotEncoder

@click.command()
@click.option("--file","-f",help="mut file from random region result")
@click.option("--outp","-o",help="outprefix of annotate file")
@click.option("--region","-r",help="functioan region type")

def main(file, region, outp):
    pwdir = os.path.dirname(os.path.abspath(__file__))
    Mut, dbdir = {}, pwdir+'/database'
    phastfile, phylopfile = dbdir+'/Conservation/hg38.phastCons100way.bw', dbdir+'/Conservation/hg38.phyloP100way.bw'
    iresfile, mirnafile, kozakfile = dbdir+'/OtherDB/IRES_literature_all-hg38.bed', dbdir+'/OtherDB/TargetScanHuman-Predicted_Target_Locations.default_predictions.hg38.bed', dbdir+'/OtherDB/Kozak_motif.txt'
    outfile = outp+'_annotation.txt'
    out = open(outfile, "w")
    if '3UTR' in region:
        out.write('Var\tSeq\tPCS\tPPS\tMFE\tmiRNA\tPAS\tAU-cls\tAU-9mer\tAU-13mer\tGU\tCU\tPumilio\tAU\tAU-elements\tGC-ratio\tAC-ratio\tAG-ratio\tAT-ratio\tTC-ratio\tTG-ratio\tA-ratio\tT-ratio\tC-ratio\tG-ratio\tAA-count\tAT-count\tAC-count\tAG-count\tTA-count\tTT-count\tTC-count\tTG-count\tCA-count\tCT-count\tCC-count\tCG-count\tGA-count\tGT-count\tGC-count\tGG-count\tA-homo\tT-homo\tC-homo\tG-homo\tAA-dimo\tAT-dimo\tAC-dimo\tAG-dimo\tTA-dimo\tTT-dimo\tTC-dimo\tTG-dimo\tCA-dimo\tCT-dimo\tCC-dimo\tCG-dimo\tGA-dimo\tGT-dimo\tGC-dimo\tGG-dimo\tsequni\n')
    elif '5UTR' in region:
        out.write('Var\tSeq\tPCS\tPPS\tMFE\tIRES\tuAUG\tuStop\tuORF\tKozak-score\tGC-ratio\tAC-ratio\tAG-ratio\tAT-ratio\tTC-ratio\tTG-ratio\tA-ratio\tT-ratio\tC-ratio\tG-ratio\tAA-count\tAT-count\tAC-count\tAG-count\tTA-count\tTT-count\tTC-count\tTG-count\tCA-count\tCT-count\tCC-count\tCG-count\tGA-count\tGT-count\tGC-count\tGG-count\tA-homo\tT-homo\tC-homo\tG-homo\tAA-dimo\tAT-dimo\tAC-dimo\tAG-dimo\tTA-dimo\tTT-dimo\tTC-dimo\tTG-dimo\tCA-dimo\tCT-dimo\tCC-dimo\tCG-dimo\tGA-dimo\tGT-dimo\tGC-dimo\tGG-dimo\tsequni\n')
    IRES, miRNA = read_ires(iresfile), read_mirna(mirnafile)
    with open(file) as f:
        for line in f:
            inf = line.strip().split('\t')
            if line.startswith('Gene'):
                continue
            gene, chrom, vpos, ref, alt, mutseq, wholeseq, strand, transcript, utrlen = inf[0], inf[1], inf[2], inf[3], inf[4], inf[5], inf[6], inf[7], inf[9], inf[10]
            if 'D' in alt or 'I' in alt or '.' in alt or 'KI' in chrom:
                continue
            var = gene+'_'+chrom+'_'+vpos+'_'+ref+'_'+alt+'_'+region
            pos = str(re.search('\(', mutseq).span()[0])
            if strand == '-':
                ref = revcomp(ref)
                alt = revcomp(alt)
            if var in Mut:
                continue
            else:
                Mut[var] = line.strip()
            flank1, flank2, vlen, vlen_alt = 10, 50, len(ref), len(alt)
            if vlen > 5 or vlen_alt > 5:
                continue
            if '3UTR' in region:
                flank1, flank2 = 10, 100
            phastcons, phylopcons = Conserv_100way(chrom, strand, vpos, phastfile, phylopfile)
            if phastcons == 'NA':
                continue
            subseq1, subseq2 = Subseq(wholeseq, flank1, vlen, pos, ref, alt), Subseq(wholeseq, flank2, vlen, pos, ref, alt)
            uAUG, uSUG, ustop, pAUG, pSUG, pstop, uORF, PAS, pPAS = Alter_num(mutseq, ref, alt, var, utrlen)
            elementStat = Element_stat(mutseq, ref, alt)
            ntStat = NT_content(subseq2, ref, alt)
            mfe = MFE_get(subseq2, ref, alt)
            oinf = ''
            if '3UTR' in region:
                mirnascore = miRNA_map(gene, chrom, vpos, miRNA)
                oinf = var+'_'+pos+'\t'+subseq2+'\t'+phastcons+'\t'+phylopcons+'\t'+mfe+'\t'+mirnascore+'\t'+PAS+'\t'+elementStat+'\t'+ntStat
            elif '5UTR' in region:
                kozak_score = Score_kozak(var, pos, subseq2, ref, alt, kozakfile)
                ires = IRES_map(gene, chrom, vpos, IRES)
                oinf = var+'_'+pos+'\t'+subseq2+'\t'+phastcons+'\t'+phylopcons+'\t'+mfe+'\t'+ires+'\t'+uAUG+'\t'+ustop+'\t'+uORF+'\t'+kozak_score+'\t'+ntStat
            out.write(oinf+'\n')
    out.close()

def IRES_map(gene, chrom, vpos, IRES):
    ires_fg = 0
    if chrom in IRES:
        for iresreg in IRES[chrom]:
            ires_st, ires_ed = iresreg.split('\t')[0], iresreg.split('\t')[1]
            if int(vpos) <= int(ires_ed) and int(vpos) >= int(ires_st):
                ires_fg += 1
    if ires_fg > 0:
        ires_fg = 1
    return str(ires_fg)

def miRNA_map(gene, chrom, vpos, miRNA):
    mirna_arry = []
    mirna_as = 0
    if gene in miRNA[chrom]:
        for mirnareg in miRNA[chrom][gene]:
            mirna_st, mirna_ed = mirnareg.split('\t')[0], mirnareg.split('\t')[1]
            if int(vpos) <= int(mirna_ed) and int(vpos) >= int(mirna_st):
                mirna_arry.append(float(miRNA[chrom][gene][mirnareg])/100)
        mirna_n, mirna_s, mirna_as = len(mirna_arry), np.median(mirna_arry), sum(mirna_arry)
    return str(mirna_as)

def read_mirna(mirnafile):
    mirna_dict = {}
    with open(mirnafile) as mif:
        for miline in mif:
            if 'chr' not in miline:
                continue
            miinf = miline.strip().split('\t')
            michr, mistart, miend, mige, miscore = miinf[0], miinf[1], miinf[2], miinf[3].split(':')[0], miinf[4]
            mirna_dict.setdefault(michr, {})
            mirna_dict[michr].setdefault(mige, {})
            mirna_dict[michr][mige][mistart+'\t'+miend] = miscore
    return mirna_dict 

def read_ires(iresfile):
    ires_dict = {}
    with open(iresfile) as iif:
        for iline in iif:
            iinf = iline.strip().split('\t')
            iichr, iistart, iiend, iitrans = iinf[0], iinf[1], iinf[2], iinf[3].split('.')[0]
            ires_dict.setdefault(iichr, {})
            ires_dict[iichr][iistart+'\t'+iiend] = 1
    return ires_dict 

def Conserv_100way(chrom, strand, pos, phastfile, phylopfile):
    phastbw, phylopbw = pyBigWig.open(phastfile), pyBigWig.open(phylopfile)
    cpos = int(pos)
    phastscore, phylopscore = 'NA', 'NA'
    if phastbw.intervals(chrom, cpos, cpos+1) and phylopbw.intervals(chrom, cpos, cpos+1):
        phastscore, phylopscore = phastbw.intervals(chrom, cpos, cpos+1)[0][2], phylopbw.intervals(chrom, cpos, cpos+1)[0][2]
    return str(phastscore), str(phylopscore)

def Subseq(wholeseq, flank, vlen, pos, ref, alt):
    u_len, d_len = math.ceil(float(flank-vlen)/2), int(float(flank-vlen)/2)
    start, end = int(pos)-u_len, int(pos)+d_len+1
    pstart, pend = 0, len(wholeseq)
    if start < pstart:
        start = pstart
    if end > pend:
        end = pend        
    preg = str(start)+'-'+str(end)
    preg1 = preg.split('-')[0]+'-'+str(int(pos))
    preg2 = str(int(pos)+len(ref))+'-'+preg.split('-')[1]
    pseq1 = wholeseq[int(preg1.split('-')[0]):int(preg1.split('-')[1])]
    pseq2 = wholeseq[int(preg2.split('-')[0]):int(preg2.split('-')[1])]
    pseq = pseq = pseq1+'('+ref+'/'+alt+')'+pseq2
    return pseq

def Alter_num(mutseq, ref, alt, var, utrlen):
    u_seq, d_seq = mutseq.split('(')[0], mutseq.split(')')[1]
    ref_seq, alt_seq = u_seq+ref+d_seq, u_seq+alt+d_seq
    flk = 2
    if '3UTR' in var:
        flk = 5
    uAUG, uCUG, uGUG, uSUG, ustop, uORF, PAS = '0', '0', '0', '0', '0', '0', '0'
    nAUG, nCUG, nGUG, nSUG, nstop, nPAS = 'NA', 'NA', 'NA', 'NA', 'NA', 'NA'
    pAUG, pCUG, pGUG, pSUG, pstop, pPAS = 'NA', 'NA', 'NA', 'NA', 'NA', 'NA'
    rAUG, rSUG, rstop, rPAS = get_element(ref_seq, flk)
    aAUG, aSUG, astop, aPAS = get_element(alt_seq, flk)
    if '5UTR' in var:
        nAUG, nSUG = str(float(len(aAUG)-len(rAUG))), str(float(len(aSUG)-len(rSUG)))
        nstop, astop_new, rstop_new = get_numstop(aAUG, rAUG, astop, rstop)
        pAUG = get_position(aAUG, rAUG)
        pSUG = get_positionv2(aSUG, rSUG, rAUG)
        pstop = get_position(astop_new, rstop_new)
        uORF = get_length(aAUG, rAUG, aSUG, rSUG, astop_new, rstop_new, utrlen)    
        if float(nAUG) > 0 and float(pAUG) == -1:
            uAUG = 1
        elif float(nAUG) < 0 and float(pAUG) == -1:
            uAUG = -1
        if float(nstop) > 0 and float(pstop) == -1:
            ustop = 1
        elif float(nstop) < 0 and float(pstop) == -1:
            ustop = -1
    elif '3UTR' in var:
        nPAS = str(float(len(aPAS)-len(rPAS)))
        pPAS = get_position(aPAS, rPAS)
        if float(nPAS) > 0 and float(pPAS) == -1:
            PAS = 1
        elif float(nPAS) < 0 and float(pPAS) == -1:
            PAS = -1
    return str(uAUG), str(uSUG), str(ustop), str(pAUG), str(pSUG), str(pstop), str(uORF), str(PAS), str(pPAS)

def get_numstop(au_ele, ru_ele, as_ele, rs_ele):
    num_s = 0
    anum_s, rnum_s = 0, 0
    as_ele_new, rs_ele_new = [], []
    if len(au_ele) > 0:
        for as_id in as_ele:
            if (as_id-au_ele[0]) > 0 and (as_id-au_ele[0]) % 3 == 0:
                anum_s += 1
                as_ele_new.append(as_id)
    if len(ru_ele) > 0:
        for rs_id in rs_ele:
            if (rs_id-ru_ele[0]) > 0 and (rs_id-ru_ele[0]) % 3 == 0:
                rnum_s += 1
                rs_ele_new.append(rs_id)
    num_s = anum_s-rnum_s
    return str(num_s), as_ele_new, rs_ele_new


def MFE_get(subseq, ref, alt):
    su_seq, sd_seq = subseq.split('(')[0], subseq.split(')')[1]
    sref_seq, salt_seq = su_seq+ref+sd_seq, su_seq+alt+sd_seq
    sref_len, salt_len = len(sref_seq), len(salt_seq)
    fafile, ofile = 'tmp.fasta', 'tmp.res'
    fa = open(fafile, "w")
    fa.write(">ref\n"+sref_seq+'\n>alt\n'+salt_seq)
    fa.close()
    os.system("RNAfold %s > %s" % (fafile, ofile))
    Mfe = {}
    with open('tmp.res') as mf:
        i = 0
        for mline in mf:
            i += 1
            if i % 3 == 1:
                #print mline.strip()
                mid = mline.strip()[1:]
            elif i % 3 == 0:
                #print mline.strip()
                mmfe = re.search("-*(\d+\.*\d*)",mline.strip()).group()
                Mfe[mid] = mmfe
    ref_mfe, alt_mfe = Mfe['ref'], Mfe['alt']
    #ref_mfe, alt_mfe = float(ref_mfe)/sref_len, float(alt_mfe)/salt_len
    #if float(ref_mfe) == 0:
    #    cmp_mfe = (float(alt_mfe)-float(ref_mfe)+1)/(float(ref_mfe)+1)
    #else:
    #    cmp_mfe = (float(alt_mfe)-float(ref_mfe))/(float(ref_mfe))
    cmp_mfe = float(alt_mfe)-float(ref_mfe)
    os.system("rm %s %s *_ss.ps" % (fafile, ofile))
    return str(cmp_mfe)

def Score_kozak(var, pos, subseq, ref, alt, kozakfile):
    OHE = convert_to_one_hot()
    usub_len, dsub_len = 7, 7
    su_seq, sd_seq = subseq.split('(')[0], subseq.split(')')[1]
    if len(su_seq) < 7:
        usub_len = len(su_seq)
    if len(sd_seq) < 7:
        dsub_len = len(sd_seq)
    ku_seq, kd_seq = su_seq[-usub_len:], sd_seq[:dsub_len]
    sref_seq, salt_seq = ku_seq+ref+kd_seq, ku_seq+alt+kd_seq
    ref_kozak, alt_kozak = Kozak_search(sref_seq), Kozak_search(salt_seq)
    tomfile, tomdir = 'tmp_vartom.txt', 'tmp_vartom_out'
    vartomfile = tomdir+'/tomtom.tsv'
    tom = open(tomfile, "w")
    tom.write("MEME version 4\n\nALPHABET= ACGT\n\nBackground letter frequencies (from uniform background):\nA 0.25000 C 0.25000 G 0.25000 T 0.25000\n")
    n = 0
    if len(ref_kozak) > 0:
        for i in range(0, len(ref_kozak)):
            n += 1
            #print ref_kozak
            tom.write("\nMOTIF "+var+'_'+pos+"_ref"+"_"+str(i+1)+" "+str(n)+"\n\nletter-probability matrix: alength= 4 w= 8 nsites= 20 E= 0\n")
            r_mtx = Mtx(ref_kozak[i], OHE)
            r_out = Print_mtx(r_mtx)
            tom.write('\n'.join(r_out))
    if len(alt_kozak) > 0:
        for j in range(0, len(alt_kozak)):
            n += 1 
            tom.write("\nMOTIF "+var+'_'+pos+"_alt"+"_"+str(j+1)+" "+str(n)+"\n\nletter-probability matrix: alength= 4 w= 8 nsites= 20 E= 0\n")
            a_mtx = Mtx(alt_kozak[j], OHE)
            a_out = Print_mtx(a_mtx)
            tom.write('\n'.join(a_out))
    tom.close()
    kzk_score = str(0)
    if len(ref_kozak) > 0 or len(alt_kozak) > 0 :
        os.system("tomtom %s %s -oc %s" % (tomfile, kozakfile, tomdir))    
        kzk_score = Similar_kozak(var+'_'+pos+"_ref", var+'_'+pos+"_alt", vartomfile)
        os.system("rm -r %s" % (tomdir))
    os.system("rm %s" % (tomfile))
    return kzk_score

def NT_content(subseq, ref, alt):
    NTinf = ['GC_ratio', 'AC_ratio', 'AG_ratio', 'AT_ratio', 'TC_ratio', 'TG_ratio', 'A_ratio', 'T_ratio', 'C_ratio', 'G_ratio', 'AA_count', 'AT_count', 'AC_count', 'AG_count', 'TA_count', 'TT_count', 'TC_count', 'TG_count', 'CA_count', 'CT_count', 'CC_count', 'CG_count', 'GA_count', 'GT_count', 'GC_count', 'GG_count', 'A_homo', 'T_homo', 'C_homo', 'G_homo', 'AA_dimo', 'AT_dimo', 'AC_dimo', 'AG_dimo', 'TA_dimo', 'TT_dimo', 'TC_dimo', 'TG_dimo', 'CA_dimo', 'CT_dimo', 'CC_dimo', 'CG_dimo', 'GA_dimo', 'GT_dimo', 'GC_dimo', 'GG_dimo', 'sequni']
    su_seq, sd_seq = subseq.split('(')[0], subseq.split(')')[1]
    sref_seq, salt_seq = su_seq+ref+sd_seq, su_seq+alt+sd_seq
    #print(subseq)
    ref_stat, alt_stat = NTstat(sref_seq), NTstat(salt_seq)
    cmpr_stat = CompareStat_nt(ref_stat, alt_stat)
    Stat_inf = Format_out(alt_stat, NTinf)
    return Stat_inf

def Element_stat(mutseq, ref, alt):
    MFname = ['AU-cls', 'AU-9mer', 'AU-13mer', 'GU', 'CU', 'Pumilio', 'AU', 'element']
    u_seq, d_seq = mutseq.split('(')[0], mutseq.split(')')[1]
    ref_seq, alt_seq = u_seq+ref+d_seq, u_seq+alt+d_seq
    ref_stat, alt_stat = MFstat(ref_seq), MFstat(alt_seq)
    cmpr_stat = CompareStat_elm(ref_stat, alt_stat)
    Stat_inf = Format_out(cmpr_stat, MFname)
    return Stat_inf

def Similar_kozak(rid, vid, tfile):
    Tom = {}
    with open(tfile) as tf:
        for tline in tf:
            if tline.startswith('Query_ID') or tline.startswith('#'):
                continue
            if len(tline.strip()) == 0:
                continue
            tinf = tline.strip().split('\t')
            tvar, mtfid, qvalue = tinf[0], tinf[1], float(tinf[5])
            if qvalue > 0.1:
                continue
            tvar_id = '_'.join(tvar.split('_')[0:-1])
            if tvar_id in Tom:
                qv = Tom[tvar_id].split('\t')[1]
                if float(qv) <= qvalue:
                    continue
            Tom[tvar_id] = mtfid+'\t'+str(qvalue)
    ref_score, alt_score = Mtf_score(rid, Tom), Mtf_score(vid, Tom)
    cmp_score = alt_score-ref_score
    return str(cmp_score)

def Mtf_score(rvid, Tom):
    score = 0
    if rvid in Tom:
        tom_mtf, tom_qvl = Tom[rvid].split('\t')[0], Tom[rvid].split('\t')[1]
        if tom_mtf == 'weak':
            score = 0.5*(1-float(tom_qvl))
        else:
            score = 1*(1-float(tom_qvl))
    return score

def get_element(sq, flk):
    Pas = ['TTTAAA','TATAAA','GATAAA','CATAAA','ATTAAA','AGTAAA','ACTAAA','AATGAA','AATATA','AATAGA','AATACA','AATAAA','AAGAAA','AACAAA','AAAAAA']
    Stop = ['TAG','TAA','TGA']
    ar_a, ar_cg, ar_s, ar_p = [], [], [], []
    for i in range(0, len(sq)-flk):
        sq_tmp = sq[i:i+flk+1]
        if sq_tmp == 'ATG':
            ar_a.append(i)
        if sq_tmp == 'CTG' or sq_tmp == 'GTG':
            ar_cg.append(i)
        if sq_tmp in Stop:
            ar_s.append(i)
        if sq_tmp in Pas:
            ar_p.append(i)
    return ar_a, ar_cg, ar_s, ar_p

def get_position(ar_ele, rr_ele):
    ps_flag = 0
    if len(ar_ele) == 0 and len(rr_ele) > 0:
        ps_flag = -1
    elif len(ar_ele) > 0 and len(rr_ele) == 0:
        ps_flag = -1
    elif len(ar_ele) > 0 and len(rr_ele) > 0:
        if ar_ele[0] != rr_ele[0]:
            ps_flag = -1
        else:
            dif1 = list(set(ar_ele)-set(rr_ele))
            dif2 = list(set(rr_ele)-set(ar_ele))
            dif = dif1 + dif2
            if len(dif) > 0:
                ps_flag = 1
    return str(ps_flag)

def get_positionv2(ar_elev2, rr_elev2, rra_elev2):
    ps_flagv2 = 0
    if len(rra_elev2) == 0:
        if len(ar_elev2) > 0 and len(rr_elev2) > 0:
            if ar_elev2[0] != rr_elev2[0]:
                ps_flagv2 = -1 
            else:
                dif1 = list(set(ar_elev2)-set(rr_elev2))
                dif2 = list(set(rr_elev2)-set(ar_elev2))
                dif = dif1 + dif2
                if len(dif) > 0:
                    ps_flagv2 = 1
        elif len(ar_elev2) + len(rr_elev2) > 0:
            ps_flagv2 = -1
    else:
        if len(ar_elev2) > 0 and len(rr_elev2) > 0:
            if ar_elev2[0] != rr_elev2[0]:
                if rra_elev2[0] <= ar_elev2[0] and rra_elev2[0] <= rr_elev2[0]:
                    ps_flagv2 = 1
                else:
                    ps_flagv2 = -1
            else:
                dif1 = list(set(ar_elev2)-set(rr_elev2))
                dif2 = list(set(rr_elev2)-set(ar_elev2))
                dif = dif1 + dif2
                if len(dif) > 0:
                    ps_flagv2 = 1
        elif len(ar_elev2) > 0 and len(rr_elev2) == 0:
            if ar_elev2[0] < rra_elev2[0]:
                ps_flagv2 = -1
            else:
                ps_flagv2 = 1
        elif len(ar_elev2) == 0 and len(rr_elev2) > 0 :
            if rr_elev2[0] < rra_elev2[0]:
                ps_flagv2 = -1
            else:
                ps_flagv2 = 1
    return str(ps_flagv2)

def get_length(a_a, r_a, a_cg, r_cg, a_s, r_s, utr_len):
    aa_len, ra_len, acg_len, rcg_len = 0, 0, 0, 0
    if len(a_a) > 0 and len(a_s) > 0:
        aa_len = a_s[0] - a_a[0]
    if len(r_a) > 0 and len(r_s) > 0:
        ra_len = r_s[0] - r_a[0]
    if len(a_a) > 0 and len(a_s) == 0:
        aa_len = float(utr_len) - a_a[0]
    if len(r_a) > 0 and len(r_s) == 0:
        ra_len = float(utr_len) - r_a[0]
    #if len(a_cg) > 0 and len(a_s) > 0:
    #    acg_len = a_s[0] - a_cg[0]
    #if len(r_s) > 0 and len(r_cg) > 0:
    #    rcg_len = r_s[0] - r_cg[0]
    #a_len_alter, cg_len_alter = float(aa_len - ra_len)/float(utr_len), float(acg_len - rcg_len)*0.5/float(utr_len)
    a_len_alter = float(aa_len - ra_len)/float(utr_len)
    len_alter = a_len_alter
    return str(len_alter)

def CompareStat_nt(R_stat, A_stat):
    C_stat = {}
    for st in R_stat:
        if A_stat[st]-R_stat[st] == 0:
            c_alter = 0
        elif R_stat[st] != 0:
            c_alter = float(A_stat[st]-R_stat[st])/(R_stat[st])
        else:
            if 'ratio' in st:
                c_alter = float(A_stat[st]-R_stat[st]+0.001)/(R_stat[st]+0.001)
            else:
                c_alter = float(A_stat[st]-R_stat[st]+0.1)/(R_stat[st]+0.1)
        #c_alter = math.log(c_alter, 2)
        C_stat[st] = c_alter
    return C_stat

def CompareStat_elm(R_stat, A_stat):
    C_stat = {}
    for st in R_stat:
        #c_alter = float(A_stat[st]-R_stat[st])/(R_stat[st]+0.1)
        #c_alter = math.log(c_alter, 2)
        c_alter = A_stat[st]-R_stat[st]
        C_stat[st] = c_alter
    return C_stat

def Format_out(S, N):
    s_inf = []
    for nc in N:
        s_inf.append(str(S[nc]))
    o_inf = '\t'.join(s_inf)
    return o_inf 

def NTstat(seq):
    NT = {'A':0, 'T':0, 'C':0, 'G':0}
    Uni = 0
    Homo = {'A':[0], 'T':[0], 'C':[0], 'G':[0]}
    Homo_len = {'A':1, 'T':1, 'C':1, 'G':1}
    Dint = {'AA':0, 'AT':0, 'AC':0, 'AG':0, 'TA':0, 'TT':0, 'TC':0, 'TG':0, 'CA':0, 'CT':0, 'CC':0, 'CG':0, 'GA':0, 'GT':0, 'GC':0, 'GG':0}
    Dimo = {'AA':[0], 'AT':[0], 'AC':[0], 'AG':[0], 'TA':[0], 'TT':[0], 'TC':[0], 'TG':[0], 'CA':[0], 'CT':[0], 'CC':[0], 'CG':[0], 'GA':[0], 'GT':[0], 'GC':[0], 'GG':[0]}
    Dimo_len = {'AA':1, 'AT':1, 'AC':1, 'AG':1, 'TA':1, 'TT':1, 'TC':1, 'TG':1, 'CA':1, 'CT':1, 'CC':1, 'CG':1, 'GA':1, 'GT':1, 'GC':1, 'GG':1}
    #print "A/T/C/G ratio start......"
    for i in range(0, len(seq)):
        if seq[i] == 'N':
            continue
        NT[seq[i]] += 1
    #print "Homo A/T/C/G start......"
    for j in range(0, len(seq)-1):
        if seq[j] == 'N':
            continue
        if seq[j+1] == seq[j]:
            Homo_len[seq[j]] += 1
        else:
            Homo[seq[j]].append(Homo_len[seq[j]])
            Homo_len[seq[j]] = 1 
        disq = seq[j]+seq[j+1]
        if 'N' in disq:
            continue
        Dint[disq] += 1
    #print "Dimo NT start......"
    for l in range(0, len(seq)-3):
        m = l
        while(m < len(seq)-3):
            if 'N' in seq[m]+seq[m+1]:
                next
            elif seq[m]+seq[m+1] == seq[m+2]+seq[m+3]:
                Dimo_len[seq[m]+seq[m+1]] += 1
            else:
                Dimo[seq[m]+seq[m+1]].append(Dimo_len[seq[m]+seq[m+1]])
                Dimo_len[seq[m]+seq[m+1]] = 1
            m += 2
    #print "Uniq start......"
    for k in range(1, len(seq)):
        if seq[k] == 'N':
            continue
        if seq[k] == seq[k-1]:
            Uni += 1

    A_homo, T_homo, C_homo, G_homo = max(Homo['A']), max(Homo['T']), max(Homo['C']), max(Homo['G'])
    AA_dimo, AT_dimo, AC_dimo, AG_dimo, TA_dimo, TT_dimo, TC_dimo, TG_dimo, CA_dimo, CT_dimo, CC_dimo, CG_dimo, GA_dimo, GT_dimo, GC_dimo, GG_dimo = max(Dimo['AA']), max(Dimo['AT']), max(Dimo['AC']), max(Dimo['AG']), max(Dimo['TA']), max(Dimo['TT']), max(Dimo['TC']), max(Dimo['TG']), max(Dimo['CA']), max(Dimo['CT']), max(Dimo['CC']), max(Dimo['CG']), max(Dimo['GA']), max(Dimo['GT']), max(Dimo['GC']), max(Dimo['GG'])
    
    Stat = {}
    Stat['A_ratio'], Stat['T_ratio'], Stat['G_ratio'], Stat['C_ratio'] = float(NT['A'])/len(seq), float(NT['T'])/len(seq), float(NT['G'])/len(seq), float(NT['C'])/len(seq)
    Stat['AT_ratio'], Stat['AG_ratio'], Stat['AC_ratio'], Stat['TG_ratio'], Stat['TC_ratio'], Stat['GC_ratio'] = float(NT['A']+NT['T'])/len(seq), float(NT['A']+NT['G'])/len(seq), float(NT['A']+NT['C'])/len(seq), float(NT['T']+NT['G'])/len(seq), float(NT['T']+NT['C'])/len(seq), float(NT['G']+NT['C'])/len(seq)
    Stat['AA_count'], Stat['AT_count'], Stat['AC_count'], Stat['AG_count'], Stat['TA_count'], Stat['TT_count'], Stat['TC_count'], Stat['TG_count'], Stat['CA_count'], Stat['CT_count'], Stat['CC_count'], Stat['CG_count'], Stat['GA_count'], Stat['GT_count'], Stat['GC_count'], Stat['GG_count'] = Dint['AA'], Dint['AT'], Dint['AC'], Dint['AG'], Dint['TA'], Dint['TT'], Dint['TC'], Dint['TG'], Dint['CA'], Dint['CT'], Dint['CC'], Dint['CG'], Dint['GA'], Dint['GT'], Dint['GC'], Dint['GG']
    Stat['A_homo'], Stat['T_homo'], Stat['C_homo'], Stat['G_homo'] = A_homo, T_homo, C_homo, G_homo
    Stat['AA_dimo'], Stat['AT_dimo'], Stat['AC_dimo'], Stat['AG_dimo'], Stat['TA_dimo'], Stat['TT_dimo'], Stat['TC_dimo'], Stat['TG_dimo'], Stat['CA_dimo'], Stat['CT_dimo'], Stat['CC_dimo'], Stat['CG_dimo'], Stat['GA_dimo'], Stat['GT_dimo'], Stat['GC_dimo'], Stat['GG_dimo'] = AA_dimo, AT_dimo, AC_dimo, AG_dimo, TA_dimo, TT_dimo, TC_dimo, TG_dimo, CA_dimo, CT_dimo, CC_dimo, CG_dimo, GA_dimo, GT_dimo, GC_dimo, GG_dimo
    Stat['sequni'] = Uni
    
    return Stat

def MFstat(seq):
    Stat = {'AU-cls':0, 'AU-9mer':0, 'AU-13mer':0, 'GU':0, 'CU':0, 'Pumilio':0, 'element':0, 'AU':0}
    for i in range(0, len(seq)):
        if i <= len(seq)-5:
            sq_au = seq[i:i+5]
            #print sq_au
            if sq_au == 'ATTTA':
                Stat['AU-cls'] += 1
                Stat['AU'] += 1
                Stat['element'] += 1
        if i <= len(seq)-9:
            sq_au_9mer = seq[i:i+9]
            id_au_9mer = AU_9mer(sq_au_9mer)
            if id_au_9mer == 1:
                Stat['AU-9mer'] += 1
                Stat['AU'] += 1
                Stat['element'] += 1
            #print "AU_9mer already"
        if i <= len(seq)-13:
            sq_au_13mer = seq[i:i+13]
            id_au_13mer = AU_13mer(sq_au_13mer)
            if id_au_13mer == 1:
                Stat['AU-13mer'] += 1
                Stat['AU'] += 1
                Stat['element'] += 1
            #print "AU_13mer already"
        if i <= len(seq)-11:
            sq_gu = seq[i:i+11]
            if sq_gu == 'TGTTTGTTTGT':
                Stat['GU'] += 1
                Stat['element'] += 1
        if i <= len(seq)-8:
            sq_pm = seq[i:i+8]
            id_pm = PM(sq_pm)
            if id_pm == 1:
                Stat['Pumilio'] += 1
                Stat['element'] += 1
            #print "Pumillio already"
        id_cu = CU(seq, i)
        if id_cu == 1:
            Stat['CU'] += 1
            Stat['element'] += 1
        #print "CU already"
        #print sq_au+'\t'+sq_au_9mer+'\t'+sq_au_13mer+'\t'+sq_gu+'\t'+sq_pm+'\t'
        #print str(i)+" base end"
    return Stat

def CU(sq_cu, ii):
    cu = 0
    if ii <= len(sq_cu)-15:
        if sq_cu[ii:ii+4] == 'CCCA' or sq_cu[ii:ii+4] == 'TCCA':
            #print "CU #1 OK start base is "+str(ii)
            #print sq_cu[ii:ii+4]
            for j in range(0, len(sq_cu)-13):
                if j > 10:
                    break
                if ii+4+j+4 > len(sq_cu):
                    break
                elif sq_cu[ii+4+j:ii+4+j+4] == 'CCCT' or sq_cu[ii+4+j:ii+4+j+4] == 'CCCA':
                    #print "CU #2 OK x is "+str(j)
                    #print sq_cu[ii+4+j:ii+4+j+4]
                    k = 0
                    while(k < len(sq_cu)-13):
                        if k > 10:
                            break
                        if ii+4+j+4+k >= len(sq_cu):
                            break
                        elif sq_cu[ii+4+j+4+k] == 'A' or sq_cu[ii+j+4+k] == 'G':
                            break
                        if ii+4+j+4+k+1+6 > len(sq_cu):
                            break
                        elif sq_cu[ii+4+j+4+k:ii+4+j+4+k+5] == 'TCCCC' or sq_cu[ii+4+j+4+k:ii+4+j+4+k+5] == 'TCTCC':
                            #print "CU #3 OK y is "+str(k)
                            #print sq_cu[ii+4+j+4+k:ii+4+j+4+k+5]
                            cu = 1
                        k += 1
    return cu

def PM(s_pm):
    pm = 0
    pm_all = ['TGTAAATA', 'TGTATATA', 'TGTAGATA', 'TGTACATA']
    if s_pm in pm_all:
        pm = 1
    return pm

def AU_9mer(s_au_9):
    au_9 = 0
    if s_au_9[0:7] == 'TTATTTA':
        w_au_9_br1 = W_ident(s_au_9[-1])
        w_au_9_br2 = W_ident(s_au_9[-2])
        if w_au_9_br1*w_au_9_br2 == 1:
            au_9 = 1
    return au_9

def AU_13mer(s_au_13):
    au_13 = 0
    if s_au_13[3:10] == 'TATTTAT':
        w_au_13_b1 = W_ident(s_au_13[0])
        w_au_13_b2 = W_ident(s_au_13[1])
        w_au_13_b3 = W_ident(s_au_13[2])
        w_au_13_br1 = W_ident(s_au_13[-1])
        w_au_13_br2 = W_ident(s_au_13[-2])
        w_au_13_br3 = W_ident(s_au_13[-3])
        if w_au_13_b1*w_au_13_b2*w_au_13_b3*w_au_13_br1*w_au_13_br2*w_au_13_br3 == 1:
            au_13 = 1
    if s_au_13[2:11] == 'ATTTATTTA':
        w_au_13_b1_2 = W_ident(s_au_13[0])
        w_au_13_b2_2 = W_ident(s_au_13[1])
        w_au_13_br1_2 = W_ident(s_au_13[-1])
        w_au_13_br2_2 = W_ident(s_au_13[-2])
        if w_au_13_b1_2*w_au_13_b2_2*w_au_13_br1_2*w_au_13_br2_2 == 1:
            au_13 = 1
    return au_13

def W_ident(nt):
    w = 0
    if nt == 'A' or nt == 'T':
        w = 1
    return w

def Print_mtx(o_mtx):
    o_str = []
    for vec in o_mtx:
        o_str.append('\t'.join(map(str,vec)))
    #print('\n'.join(o_str))
    return o_str
            

def Mtx(kseq, ohe):
    sq_array = []
    for k in kseq:
        sq_array.append(list(k))
    mtf_array = ohe.transform(sq_array).toarray()
    return mtf_array

def convert_to_one_hot():
    OHE = OneHotEncoder(handle_unknown='ignore')
    OHE.fit([['A'],['C'],['G'],['T']])
    return OHE

def Kozak_search(mseq):
    AUG_motif =  []
    aug_pos = []
    if 'ATG' in mseq:
        for match in re.finditer('ATG',mseq):
            aug_pos.append(match.start())
    for ap in aug_pos:
        if int(ap) >= 3 and len(mseq) >= int(ap)+5:
            ug_motif = mseq[int(ap)-3:int(ap)]+'ATG'+mseq[int(ap)+3:int(ap)+5]
            AUG_motif.append(ug_motif)
    return AUG_motif

def complement(seq2):
    trantab = seq2.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh')
    return seq2.translate(trantab)

def revcomp(seq1):
    return complement(seq1)[::-1]

if __name__ == "__main__":
    main()
