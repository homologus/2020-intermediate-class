#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Couples overlapping sequences of text (and DNA)
# back together after sequencing. Programmed by Samuel Davenport

from Bio.Seq import Seq
from Bio import SeqIO
from random import randint, choice
import os
from threading import Thread

NAME = 'Recoupler'
__version__ = '1.0.0'

DATA = 'possible-copyin, tely-suggests-a, our-notice-that, ice-that-the-sp, ible-copying-me,uggests-a-possi, -suggests-a-pos, anism-for-the-g, lated-immediate, ped-our-notice-,ly-suggests-a-p, e-copying-mecha, ying-mechanism-, ts-a-possible-c, ossible-copying,pying-mechanism, the-genetic-mat, -immediately-su, pecific-pairing, uggests-a-possi,e-copying-mecha, copying-mechani, tely-suggests-a, aped-our-notice, immediately-sug,ated-immediatel, ng-we-have-post, ing-mechanism-f, mediately-sugge, iately-suggests,suggests-a-poss, ur-notice-that-, -have-postulate, ur-notice-that-, fic-pairing-we-,-we-have-postul, e-genetic-mater, ostulated-immed, lated-immediate, e-genetic-mater,ostulated-immed, ed-immediately-, ce-that-the-spe, notice-that-the, -pairing-we-hav,ing-we-have-pos, -for-the-geneti, as-not-escaped-, t-has-not-escap, -postulated-imm,e-postulated-im, t-the-specific-, ests-a-possible, -not-escaped-ou, immediately-sug,e-that-the-spec, ossible-copying, iately-suggests, otice-that-the-, ng-we-have-post,-pairing-we-hav, has-not-escaped, e-have-postulat, ely-suggests-a-, ely-suggests-a-,e-specific-pair, ble-copying-mec, -genetic-materi, -our-notice-tha, mmediately-sugg,pying-mechanism, -our-notice-tha, notice-that-the, caped-our-notic, tice-that-the-s,-that-the-speci, e-copying-mecha, -copying-mechan, possible-copyin, ve-postulated-i,ted-immediately, -postulated-imm, ng-mechanism-fo, specific-pairin, e-specific-pair,ice-that-the-sp, -for-the-geneti, ulated-immediat, -the-specific-p, escaped-our-not,ediately-sugges, at-the-specific, ped-our-notice-, -that-the-speci, ests-a-possible,ggests-a-possib, ecific-pairing-, not-escaped-our, sm-for-the-gene, ve-postulated-i,for-the-genetic, the-genetic-mat, ely-suggests-a-, ely-suggests-a-, -copying-mechan,our-notice-that, ing-we-have-pos, ave-postulated-, stulated-immedi, -not-escaped-ou,caped-our-notic, e-that-the-spec, a-possible-copy, notice-that-the, he-specific-pai,e-have-postulat, ic-pairing-we-h, -possible-copyi, ng-we-have-post, It-has-not-esca,r-the-genetic-m, the-genetic-mat, stulated-immedi, g-mechanism-for, sible-copying-m,that-the-specif, e-specific-pair, has-not-escaped, ve-postulated-i, ic-pairing-we-h,ur-notice-that-, ing-we-have-pos, -notice-that-th, ible-copying-me, otice-that-the-,ve-postulated-i, -postulated-imm, iately-suggests, d-our-notice-th, cific-pairing-w,diately-suggest, pairing-we-have, -postulated-imm, t-the-specific-, g-mechanism-for,stulated-immedi, ur-notice-that-, a-possible-copy, e-postulated-im, anism-for-the-g,or-the-genetic-, anism-for-the-g, we-have-postula, ely-suggests-a-, r-notice-that-t,t-escaped-our-n, -copying-mechan, chanism-for-the, t-escaped-our-n, we-have-postula,specific-pairin, r-notice-that-t, hanism-for-the-, ce-that-the-spe, copying-mechani,pying-mechanism, for-the-genetic, ecific-pairing-, ostulated-immed, ly-suggests-a-p,t-has-not-escap, we-have-postula, otice-that-the-, ng-we-have-post, has-not-escaped,ice-that-the-sp, ble-copying-mec, tice-that-the-s, r-the-genetic-m, as-not-escaped-,the-genetic-mat, e-specific-pair, immediately-sug, aped-our-notice, mediately-sugge,ice-that-the-sp, has-not-escaped, c-pairing-we-ha, -the-genetic-ma, diately-suggest,that-the-specif, mmediately-sugg, -copying-mechan, ying-mechanism-, pying-mechanism,-for-the-geneti, e-postulated-im, g-we-have-postu, y-suggests-a-po, fic-pairing-we-,hat-the-specifi, scaped-our-noti, pying-mechanism, specific-pairin, pying-mechanism,as-not-escaped-, tulated-immedia, not-escaped-our, -copying-mechan, lated-immediate,ce-that-the-spe, ism-for-the-gen, ulated-immediat, ng-we-have-post, chanism-for-the,anism-for-the-g, ic-pairing-we-h, -genetic-materi, ble-copying-mec, diately-suggest,possible-copyin, y-suggests-a-po, -possible-copyi, possible-copyin, pairing-we-have,m-for-the-genet, postulated-imme, copying-mechani, ice-that-the-sp, sible-copying-m,possible-copyin, ately-suggests-, ssible-copying-, opying-mechanis, escaped-our-not'
##DATA = 'GRTGTAYLMNKVVNV, SQHRSQKNVSFITYG, SPSSSYKSINDPLFH, VCDILCHWCKRNVGW, DDQQYKEGKFILELK,TYGCRHCKTHLSSSF, FITYGCRHCKTHLSS, QYKEGKFILELKNIC, SSNDDQQYKEGKFIL, IISRDYRGRTGTAYL,LSSSFQIISRDYRGR, SSPSSSYKSINDPLF, KYLQSSNDDQQYKEG, YKEGKFILELKNICK, TYGCRHCKTHLSSSF,YKEGKFILELKNICK, SSSFQIISRDYRGRT, LVCDILCHWCKRNVG, VNVVEGKVEQRRMLT, SSSYKSINDPLFHSQ,SYKSINDPLFHSQHR, HLSSSFQIISRDYRG, LMNKVVNVVEGKVEQ, RDYRGRTGTAYLMNK, KRNVGWKYLQSSNDD,YSIYIENPLSSPSSS, RHCKTHLSSSFQIIS, GDYLVCDILCHWCKR, SSFQIISRDYRGRTG, QYKEGKFILELKNIC,SSSYKSINDPLFHSQ, VVNVVEGKVEQRRML, THLSSSFQIISRDYR, RMLTGDYLVCDILCH, PLSSPSSSYKSINDP,SSFQIISRDYRGRTG, THLSSSFQIISRDYR, YLMNKVVNVVEGKVE, ITYGCRHCKTHLSSS, SFITYGCRHCKTHLS,VVNVVEGKVEQRRML, SSNDDQQYKEGKFIL, TGTAYLMNKVVNVVE, ITYGCRHCKTHLSSS, DYLVCDILCHWCKRN,QIISRDYRGRTGTAY, DYLVCDILCHWCKRN, TGTAYLMNKVVNVVE, DQQYKEGKFILELKN, GKVEQRRMLTGDYLV,WKYLQSSNDDQQYKE, SSNDDQQYKEGKFIL, DILCHWCKRNVGWKY, TYGCRHCKTHLSSSF, MLTGDYLVCDILCHW,DDQQYKEGKFILELK, YLVCDILCHWCKRNV, YLQSSNDDQQYKEGK, EGKFILELKNICKCT, PLSSPSSSYKSINDP,YLVCDILCHWCKRNV, RMLTGDYLVCDILCH, QYKEGKFILELKNIC, CRHCKTHLSSSFQII, NKVVNVVEGKVEQRR,EQRRMLTGDYLVCDI, SSSFQIISRDYRGRT, YRGRTGTAYLMNKVV, SRDYRGRTGTAYLMN, SPSSSYKSINDPLFH,GKVEQRRMLTGDYLV, VSFITYGCRHCKTHL, SIYIENPLSSPSSSY, KSINDPLFHSQHRSQ, NPLSSPSSSYKSIND,SQKNVSFITYGCRHC, GTAYLMNKVVNVVEG, FHSQHRSQKNVSFIT, AYLMNKVVNVVEGKV, YKSINDPLFHSQHRS,FHSQHRSQKNVSFIT, RRMLTGDYLVCDILC, RMLTGDYLVCDILCH, TAYLMNKVVNVVEGK, CHWCKRNVGWKYLQS,TAYLMNKVVNVVEGK, KRNVGWKYLQSSNDD, SYKSINDPLFHSQHR, KVEQRRMLTGDYLVC, SINDPLFHSQHRSQK,SQHRSQKNVSFITYG, EGKFILELKNICKCT, MLTGDYLVCDILCHW, VGWKYLQSSNDDQQY, SSPSSSYKSINDPLF,CKRNVGWKYLQSSND, RGRTGTAYLMNKVVN, LMNKVVNVVEGKVEQ, FQIISRDYRGRTGTA, LFHSQHRSQKNVSFI,SSSFQIISRDYRGRT, AYLMNKVVNVVEGKV, DYRGRTGTAYLMNKV, MNKVVNVVEGKVEQR, MGLRYSIYIENPLSS,KSINDPLFHSQHRSQ, LSSPSSSYKSINDPL, NDPLFHSQHRSQKNV, PSSSYKSINDPLFHS, NDDQQYKEGKFILEL,CKRNVGWKYLQSSND, SSFQIISRDYRGRTG, VCDILCHWCKRNVGW, RDYRGRTGTAYLMNK, QKNVSFITYGCRHCK,LSSPSSSYKSINDPL, LSSSFQIISRDYRGR, KVVNVVEGKVEQRRM, QKNVSFITYGCRHCK, YGCRHCKTHLSSSFQ,FITYGCRHCKTHLSS, HCKTHLSSSFQIISR, WKYLQSSNDDQQYKE, CDILCHWCKRNVGWK, RTGTAYLMNKVVNVV,SPSSSYKSINDPLFH, SSSFQIISRDYRGRT, SYKSINDPLFHSQHR, TAYLMNKVVNVVEGK, ENPLSSPSSSYKSIN,GLRYSIYIENPLSSP, SIYIENPLSSPSSSY, DILCHWCKRNVGWKY, KSINDPLFHSQHRSQ, KEGKFILELKNICKC,LFHSQHRSQKNVSFI, YKEGKFILELKNICK, KEGKFILELKNICKC, LSSPSSSYKSINDPL, YLQSSNDDQQYKEGK,SSSYKSINDPLFHSQ, NDDQQYKEGKFILEL, LMNKVVNVVEGKVEQ, VNVVEGKVEQRRMLT, ISRDYRGRTGTAYLM,RDYRGRTGTAYLMNK, LSSSFQIISRDYRGR, SSPSSSYKSINDPLF, NVVEGKVEQRRMLTG, VGWKYLQSSNDDQQY,YGCRHCKTHLSSSFQ, YLVCDILCHWCKRNV, YKEGKFILELKNICK, LQSSNDDQQYKEGKF, GDYLVCDILCHWCKR,RDYRGRTGTAYLMNK, VVEGKVEQRRMLTGD, YKSINDPLFHSQHRS, KVVNVVEGKVEQRRM, YSIYIENPLSSPSSS,WKYLQSSNDDQQYKE, QHRSQKNVSFITYGC, QHRSQKNVSFITYGC, YLQSSNDDQQYKEGK, SSYKSINDPLFHSQH,SFITYGCRHCKTHLS, VCDILCHWCKRNVGW, LTGDYLVCDILCHWC, TGTAYLMNKVVNVVE, RDYRGRTGTAYLMNK,SSNDDQQYKEGKFIL, EQRRMLTGDYLVCDI, RYSIYIENPLSSPSS, DILCHWCKRNVGWKY, YGCRHCKTHLSSSFQ,IYIENPLSSPSSSYK, RSQKNVSFITYGCRH, DQQYKEGKFILELKN, TGTAYLMNKVVNVVE, QHRSQKNVSFITYGC,SSYKSINDPLFHSQH, NPLSSPSSSYKSIND, YKSINDPLFHSQHRS, EQRRMLTGDYLVCDI, DQQYKEGKFILELKN'
##DATA = 'TGYKLKNGRSCKSCW, LPWTYPPRFYCSKCG, RFYCSKCGNTGYKLK, YTMPVYTNAWQGNRP, PLYVQPGDPRLGGVL,YYTNYTMPVYTNAWQ, GVLCGECRGSGRTRF, AWQGNRPLYVQPGDP, LYVQPGDPRLGGVLC, LPWTYPPRFYCSKCG,QGNRPLYVQPGDPRL, TPTSSQPRPPPRPQQ, PHQRPSTMPATSSSQ, NLAQGHQSRPHQRPS, QGHQSRPHQRPSTMP,APQNNVVSAPTYYTN, YYTNYTMPVYTNAWQ, TYAHSHSYTPTSSQP, PRFYCSKCGNTGYKL, WTYPPRFYCSKCGNT,QNNVVSAPTYYTNYT, TSSSQTYAHSHSYTP, KLKNGRSCKSCWRRF, QPPRPPRPAANLAQG, CSKCGNTGYKLKNGR,EEERLQSQPPRPPRP, AANLAQGHQSRPHQR, PATSSSQTYAHSHSY, QTYAHSHSYTPTSSQ, TPTSSQPRPPPRPQQ,ERLQSQPPRPPRPAA, GGVLCGECRGSGRTR, TRFLLDEDICPLCHG, QNPSLPWTYPPRFYC, QGHQSRPHQRPSTMP,PPRPAANLAQGHQSR, RSCKSCWRRFAPQNN, GNTGYKLKNGRSCKS, KEEERLQSQPPRPPR, DTHDDELPSYEDVIK,PTYYTNYTMPVYTNA, NRPLYVQPGDPRLGG, TMPVYTNAWQGNRPL, NTGYKLKNGRSCKSC, FLLDEDICPLCHGVG,DTHDDELPSYEDVIK, WRRFAPQNNVVSAPT, GNRPLYVQPGDPRLG, RLQSQPPRPPRPAAN, PHQRPSTMPATSSSQ,NPSLPWTYPPRFYCS, PRPQQNPSLPWTYPP, APQNNVVSAPTYYTN, RPQQNPSLPWTYPPR, SKCGNTGYKLKNGRS,YTPTSSQPRPPPRPQ, RGSGRTRFLLDEDIC, FLLDEDICPLCHGVG, ERLQSQPPRPPRPAA, EEERLQSQPPRPPRP,CGECRGSGRTRFLLD, TPTSSQPRPPPRPQQ, NRPLYVQPGDPRLGG, WRRFAPQNNVVSAPT, HSHSYTPTSSQPRPP,APTYYTNYTMPVYTN, PSYEDVIKEEERLQS, RPQQNPSLPWTYPPR, SRPHQRPSTMPATSS, GGVLCGECRGSGRTR,YPPRFYCSKCGNTGY, QRPSTMPATSSSQTY, LGGVLCGECRGSGRT, ATSSSQTYAHSHSYT, LYVQPGDPRLGGVLC,QNNVVSAPTYYTNYT, TYPPRFYCSKCGNTG, SSQPRPPPRPQQNPS, PHQRPSTMPATSSSQ, KCGNTGYKLKNGRSC,QPRPPPRPQQNPSLP, ATSSSQTYAHSHSYT, LPWTYPPRFYCSKCG, FYCSKCGNTGYKLKN, TNAWQGNRPLYVQPG,CGECRGSGRTRFLLD, LCHGVGRIITQPQRY, YTMPVYTNAWQGNRP, TMPATSSSQTYAHSH, AANLAQGHQSRPHQR,VVSAPTYYTNYTMPV, TGYKLKNGRSCKSCW, NPSLPWTYPPRFYCS, PGDPRLGGVLCGECR, GRTRFLLDEDICPLC,RPLYVQPGDPRLGGV, PVYTNAWQGNRPLYV, YEDVIKEEERLQSQP, NYTMPVYTNAWQGNR, SKCGNTGYKLKNGRS,KSCWRRFAPQNNVVS, MPATSSSQTYAHSHS, KEEERLQSQPPRPPR, SQPRPPPRPQQNPSL, LAQGHQSRPHQRPST,ECRGSGRTRFLLDED, SQPRPPPRPQQNPSL, CRGSGRTRFLLDEDI, YCSKCGNTGYKLKNG, ECRGSGRTRFLLDED,DELPSYEDVIKEEER, APQNNVVSAPTYYTN, PSTMPATSSSQTYAH, QPRPPPRPQQNPSLP, LPWTYPPRFYCSKCG,LGGVLCGECRGSGRT, PRPPRPAANLAQGHQ, AANLAQGHQSRPHQR, STMPATSSSQTYAHS, PRPQQNPSLPWTYPP,YTNYTMPVYTNAWQG, CPLCHGVGRIITQPQ, QSRPHQRPSTMPATS, LQSQPPRPPRPAANL, RFAPQNNVVSAPTYY,DPRLGGVLCGECRGS, QSRPHQRPSTMPATS, HDDELPSYEDVIKEE, YEDVIKEEERLQSQP, TNYTMPVYTNAWQGN,PSTMPATSSSQTYAH, YPPRFYCSKCGNTGY, ERLQSQPPRPPRPAA, NRPLYVQPGDPRLGG, TYAHSHSYTPTSSQP,RFAPQNNVVSAPTYY, LDEDICPLCHGVGRI, CKSCWRRFAPQNNVV, QPPRPPRPAANLAQG, FAPQNNVVSAPTYYT,LLDEDICPLCHGVGR, QNNVVSAPTYYTNYT, VYTNAWQGNRPLYVQ, DELPSYEDVIKEEER, QQNPSLPWTYPPRFY,PATSSSQTYAHSHSY, ICPLCHGVGRIITQP, TMPVYTNAWQGNRPL, VYTNAWQGNRPLYVQ, PLYVQPGDPRLGGVL,SQTYAHSHSYTPTSS, NRPLYVQPGDPRLGG, TSSQPRPPPRPQQNP, EDICPLCHGVGRIIT, HDDELPSYEDVIKEE,RRFAPQNNVVSAPTY, GHQSRPHQRPSTMPA, DDELPSYEDVIKEEE, PTSSQPRPPPRPQQN, DICPLCHGVGRIITQ,PTSSQPRPPPRPQQN, RRFAPQNNVVSAPTY, RRFAPQNNVVSAPTY, TNYTMPVYTNAWQGN, KLKNGRSCKSCWRRF,GNTGYKLKNGRSCKS, PPPRPQQNPSLPWTY, QPPRPPRPAANLAQG, QSRPHQRPSTMPATS, NVVSAPTYYTNYTMP,VQPGDPRLGGVLCGE, KNGRSCKSCWRRFAP, CGNTGYKLKNGRSCK, PAANLAQGHQSRPHQ, TPTSSQPRPPPRPQQ,QSRPHQRPSTMPATS, SSQTYAHSHSYTPTS, YTPTSSQPRPPPRPQ, GNTGYKLKNGRSCKS, YTNAWQGNRPLYVQP,CGECRGSGRTRFLLD, VVSAPTYYTNYTMPV, GSGRTRFLLDEDICP, KDTHDDELPSYEDVI, YVQPGDPRLGGVLCG,GGVLCGECRGSGRTR, DEDICPLCHGVGRII, SCKSCWRRFAPQNNV, ATSSSQTYAHSHSYT, RPPRPAANLAQGHQS,HSHSYTPTSSQPRPP, FYCSKCGNTGYKLKN, PSLPWTYPPRFYCSK, LQSQPPRPPRPAANL, PRPPRPAANLAQGHQ,PQQNPSLPWTYPPRF, PHQRPSTMPATSSSQ, LPWTYPPRFYCSKCG, HQSRPHQRPSTMPAT, GVLCGECRGSGRTRF,TMPATSSSQTYAHSH, PRPPPRPQQNPSLPW, EERLQSQPPRPPRPA, PATSSSQTYAHSHSY, CGECRGSGRTRFLLD,QNNVVSAPTYYTNYT, APQNNVVSAPTYYTN, SQPRPPPRPQQNPSL, RPPRPAANLAQGHQS, QTYAHSHSYTPTSSQ,YTMPVYTNAWQGNRP, PGDPRLGGVLCGECR, FYCSKCGNTGYKLKN, SKDTHDDELPSYEDV, PVYTNAWQGNRPLYV,RPHQRPSTMPATSSS, YTNYTMPVYTNAWQG, LLDEDICPLCHGVGR, QPRPPPRPQQNPSLP, KCGNTGYKLKNGRSC,RLGGVLCGECRGSGR, HQSRPHQRPSTMPAT, SYTPTSSQPRPPPRP, GRTRFLLDEDICPLC, SQPPRPPRPAANLAQ,SHSYTPTSSQPRPPP, DELPSYEDVIKEEER, KNGRSCKSCWRRFAP, NYTMPVYTNAWQGNR, APTYYTNYTMPVYTN,EERLQSQPPRPPRPA, GGVLCGECRGSGRTR, PWTYPPRFYCSKCGN, FYCSKCGNTGYKLKN, PVYTNAWQGNRPLYV,EDVIKEEERLQSQPP, THDDELPSYEDVIKE, VIKEEERLQSQPPRP, CSKCGNTGYKLKNGR, DICPLCHGVGRIITQ,LKNGRSCKSCWRRFA, QSQPPRPPRPAANLA, RPHQRPSTMPATSSS, SCKSCWRRFAPQNNV, RLQSQPPRPPRPAAN,PPRPAANLAQGHQSR, PRPPRPAANLAQGHQ, YEDVIKEEERLQSQP, PRPPPRPQQNPSLPW, WRRFAPQNNVVSAPT,LAQGHQSRPHQRPST, SHSYTPTSSQPRPPP, TNYTMPVYTNAWQGN, THDDELPSYEDVIKE, SKDTHDDELPSYEDV,PPRPAANLAQGHQSR, QPRPPPRPQQNPSLP, RFYCSKCGNTGYKLK, GSGRTRFLLDEDICP, APQNNVVSAPTYYTN,RGSGRTRFLLDEDIC, PVYTNAWQGNRPLYV, RGSGRTRFLLDEDIC, PQNNVVSAPTYYTNY, YTPTSSQPRPPPRPQ,PAANLAQGHQSRPHQ, SQTYAHSHSYTPTSS, YYTNYTMPVYTNAWQ, PHQRPSTMPATSSSQ, AHSHSYTPTSSQPRP,QTYAHSHSYTPTSSQ, VSAPTYYTNYTMPVY, CKSCWRRFAPQNNVV, RLGGVLCGECRGSGR, LGGVLCGECRGSGRT,YVQPGDPRLGGVLCG, RPLYVQPGDPRLGGV, SHSYTPTSSQPRPPP, SQPRPPPRPQQNPSL, YYTNYTMPVYTNAWQ,NPSLPWTYPPRFYCS, GRTRFLLDEDICPLC, YTNAWQGNRPLYVQP, SSSQTYAHSHSYTPT, YTNYTMPVYTNAWQG,PTYYTNYTMPVYTNA, TYAHSHSYTPTSSQP, SSQPRPPPRPQQNPS, PRLGGVLCGECRGSG, RFYCSKCGNTGYKLK,KEEERLQSQPPRPPR, TPTSSQPRPPPRPQQ, QQNPSLPWTYPPRFY, PRPAANLAQGHQSRP, SYTPTSSQPRPPPRP,ERLQSQPPRPPRPAA, QPGDPRLGGVLCGEC, SKDTHDDELPSYEDV, SQPRPPPRPQQNPSL, SLPWTYPPRFYCSKC,SCKSCWRRFAPQNNV, WRRFAPQNNVVSAPT, SYEDVIKEEERLQSQ, GRTRFLLDEDICPLC, QSRPHQRPSTMPATS,HSHSYTPTSSQPRPP, HDDELPSYEDVIKEE, PVYTNAWQGNRPLYV, LPWTYPPRFYCSKCG, YVQPGDPRLGGVLCG,CKSCWRRFAPQNNVV, GSGRTRFLLDEDICP, GDPRLGGVLCGECRG, PRLGGVLCGECRGSG, GECRGSGRTRFLLDE,PPRFYCSKCGNTGYK, HQSRPHQRPSTMPAT, VYTNAWQGNRPLYVQ, LAQGHQSRPHQRPST, NPSLPWTYPPRFYCS,RPAANLAQGHQSRPH, TGYKLKNGRSCKSCW, SQPPRPPRPAANLAQ, LLDEDICPLCHGVGR, EEERLQSQPPRPPRP'
DATA = [i[1:] if not ' ' in i[1:] else i[1:][:i.index(' ')] for i in DATA.split(',')]
temp = []
[temp.append(i) for i in DATA if not i in temp]
DATA = list(temp)

def get_frames(dna_seq, n=3):
    """Returns proteins with a length at least <l> long from frame shifts <n>"""
    # For each frame shift,
    data = []
    for shift in range(abs(n) % 4):
        leng = len(dna_seq[shift:])
        mod = leng % 3
        #print(leng, mod)
        # Get the frame shift of the DNA sequence and
        # translate data to proteins
        proteins_data = dna_seq[shift:-mod].translate()
        # Split the proteins by their end codons
        proteins = str(proteins_data).split('*')
        # For each protein we found,
##        for protein in proteins:
##            # If the protein has at least 100 amino acids,
##            if len(protein) >= l:
##                # Write the protein to the file
##                file.write(protein+'\n')
        data += proteins
    return [protein for protein in data if len(protein) >= 5]

def gen_random_proper_seq(length, a, t, g, c):#psat=0.5, psgc=0.5):
    if sum([a, t, g, c]) != 1:
        raise ArithmeticError('Sum of perentages of A, T, G, and C is not equal to 100 percent!')
    a = ['A' for i in range(round(length * a))]
    t = ['T' for i in range(round(length * t))]
    g = ['G' for i in range(round(length * g))]
    c = ['C' for i in range(round(length * c))]
    unrand = sum([a, t, g, c], [])
    # Help free up memory
    del a, t, g, c
##    at = length * psat
##    gc = length * psgc
##    unrand = sum([['A', 'T'] for i in range(round(at/2))], []) + sum([['G', 'C'] for i in range(round(gc/2))], [])
    seq = []
    for i in range(length):
        rng = (0, len(unrand)-1)
        idx = randint(*rng)
        seq.append(unrand[idx])
        del unrand[idx]
    #print(Seq(''.join(seq)))
    return Seq(''.join(seq))

def recouple(lst, scanLen=8, start=''):
    """Program that will re-couple text from a list that has been broken into many overlapping parts."""
    # Get all the seperate chuncks of data in a list, no repeats.
    data = []
    [data.append(i) for i in lst if not i in data]
    if scanLen < 4:
        scanLen = min(round(min([len(i) for i in data])/2), 4)
        print('scanLen =',scanLen)
    # Get all the words into a dictionary
    wordsLst = {}
    beginsings = {}
    # For each thing to test in the dictionary
    for test in data:
        # If we haven't seen this word before, set it as a possible
        # beginning value
        if not test[0:scanLen] in beginsings.keys():
            beginsings[test[0:scanLen]] = True
        # For each index point our test string,
        for i in range(len(test)-scanLen):
            # Get the word it could be
            word = test[i:i+scanLen]
            # Add the word to the words dict, pointing to the last letter of
            # our test string.
            wordsLst[word] = test[i+scanLen]
            # If it's not the beginning
            if i > 0:
                beginsings[word] = False
    # If the starting word was not defined,
    if not start:
        # Set it to a random valid word
        for word in beginsings.keys():
            if beginsings[word]:
                test = str(word)
                break
    else:
        # Otherwise, try to use it
        test = str(start)
        cstart = False
        yes = False
        # If it's invalid, break
        if not test in wordsLst.keys():
            for word in wordsLst.keys():
                if test in word:
                    test = word
                    yes = True
                    break
        if not yes:
            raise ValueError('Start word argument is invalid, does not exist in word list with scanLen. Try having start word len be equal to scanLen argument number.')
    # Free up some memory
    del beginsings
    # Get data that points to eachother
    string = test
    count = len(wordsLst.keys())*scanLen
    while count:
        # If the test word is valid
        if test in wordsLst.keys():
            # Add it to the string
            string += wordsLst[test]
            test = test[1:scanLen]+wordsLst[test]
            count -= 1
        else:# Otherwise, decrement scan by one
            count = 0
##    # Once we have data that points to other bits, put them together properly
##    data = str(string[0])
##    # For each bit of data to add,
##    for i in string[1:]:
##        # For all possible index positions in the new bit of data,
##        for ii in range(len(i)):
##            # Get the current data's tail with ii
##            piece = data[-ii:]
##            # If the tail of this index is the start of this bit of data,
##            if piece == i[:len(piece)]:
##                # Combine the two chucks properly and stop checking indexes for this new chunk.
##                data = data[:len(data)-ii] + i
##                break
    return string

def run():
    print('Running...')
    print('"',recouple(DATA),'"')
##    data = {len(i):i for i in [recouple(DATA, DATA[i][:6], 6) for i in range(len(DATA))]}
##    print(max(data.keys()), data[max(data.keys())])
    print('Done')

if __name__ == '__main__':
    run()
