#!/usr/bin/env python
from __future__ import division

import sys
#import gzip
import difflib
import argparse

# import regex
import pandas as pd

def is_novel_peptide(peptide, assigned_peptides, epitope_len):
    for assigned_peptide in assigned_peptides:        
        matcher = difflib.SequenceMatcher(None, peptide, assigned_peptide)
        match = matcher.find_longest_match(0, len(peptide), 0, len(assigned_peptide))
        if match.size >= epitope_len:
            return False
    return True

# def is_novel_peptide(peptide, assigned_peptides, epitope_len, max_mm):
#     for i in range(len(peptide) - epitope_len):
#         epitope = peptide[i:i + epitope_len]
#         fmt = r'({epitope}){{s<={max_mm}}}'.format(epitope=epitope, max_mm=max_mm)
#         r = regex.compile(fmt)
#         for assigned_peptide in assigned_peptides:
#             if r.search(assigned_peptide) is not None:
#                 return False
#     return True


def calc_virus_scores(series, level, epitope_len):
    # order viruses by decreasing number of hits
    nhits_per_virus = series.groupby(level=level).sum()

    #print(nhits_per_virus.head())
    #Species
    #Absidia glauca (Pin mould)              0
    #Acidocella sp. MX-AZ02                  0
    #Adeno-associated dependoparvovirus A    0
    #Adeno-associated virus                  0
    #Adeno-associated virus VR-355           0
    #Name: S148, dtype: int64

    # initialize all scores to 0
    virus_scores = pd.Series(index=nhits_per_virus.index, name=series.name).fillna(0).astype(int)

    #print(virus_scores.head())
    #Species
    #Absidia glauca (Pin mould)              0
    #Acidocella sp. MX-AZ02                  0
    #Adeno-associated dependoparvovirus A    0
    #Adeno-associated virus                  0
    #Adeno-associated virus VR-355           0
    #Name: S148, dtype: int64

    assigned_peptides = set()
    grouped = series[series].groupby(level=level)
    #for virus in nhits_per_virus[nhits_per_virus > 0].order(ascending=False).index:

    #print(nhits_per_virus[nhits_per_virus > 0].sort_values(ascending=False).head())
    #Species
    #Influenza A virus           100
    #Human herpesvirus 4          84
    #Enterovirus B                44
    #Rhinovirus A                 24
    #Streptococcus pneumoniae     23
    #Name: S148, dtype: int64

    for virus in nhits_per_virus[nhits_per_virus > 0].sort_values(ascending=False).index:
        virus_hits = grouped.get_group(virus)

        score = 0
        peptides = virus_hits.index.get_level_values('peptide')

        #print(peptides[0:2])
        #Index(['PGDTIIFEANGNLIAPWYAFALSRGFGSGIITSNASMGECDAKCQTPQGAINSSLP', 
        #  'QIASNENMETIDSITLELRSKYWAIRTRSGGNTNKQRASAGQISVQPTFSVQRNLP'], dtype='object', name='peptide')

        for peptide in peptides:
            if is_novel_peptide(peptide, assigned_peptides, epitope_len):
            # if is_novel_peptide(peptide, assigned_peptides, epitope_len, max_mm):
                score += 1
                assigned_peptides.add(peptide)
        #print(score)
        #10
        #42
        #5
        #2
        #11
        #print(assigned_peptides)
        #{'QDLPGKSNSTATLCLGHHAVPNGTLVKTITDDQIEVTNATKLVQSSSTGRICNNPH', 'RSSSSYCRNPNNEKGTHGVKGWAFDDGNDVWMGRTISEDSRSGYETFKVIGGWSTP', 
        # 'GPDNGAVAVLKYNGIITDTIKSWRNNIMRTQESECACVNGSCFTVMTDGPSNGQAS', 'NGMVTRFESLKIYRDSLGETVMRMGDLHYLQSRNEKWRDQLGQKFEEIRWLIEEMR', 
        # 'NSQVYSLIRPNENPAHKSQLVWMACHSAAFEDLRLLSFIRGTKVCPRGKLSTRGVQ', 'QIASNENMETIDSITLELRSKYWAIRTRSGGNTNKQRASAGQISVQPTFSVQRNLP', 
        # 'SNLNDTTYQRTRALVRTGMDPRMCSLMQGSTLPRRSGAAGAAVKGVGTMVLELIRM', 'RAVGKCPRYVKQKSLLLATGMKNVPEIPKKREKRGLFGAIAGFIENGWEGLVDGWY', 
        # 'PGDTIIFEANGNLIAPWYAFALSRGFGSGIITSNASMGECDAKCQTPQGAINSSLP', 'SILNTSQRGILEDGQMYQKCCNLFEKFFPSSSYRRPVGISSMVEAMVSRARIDARI'}

        virus_scores[virus] = score
    return virus_scores


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("hits")
    parser.add_argument("oligo_metadata")
    #parser.add_argument("beads_nhits")
    #parser.add_argument("samps_nhits")
    parser.add_argument("level")
    parser.add_argument("epitope_len", type=int)
    # parser.add_argument("max_mismatch", type=int)
    args = parser.parse_args()

    hits = pd.read_csv(args.hits, index_col=0)

    lib = pd.read_csv(args.oligo_metadata,low_memory=False)
    # DtypeWarning: Columns (7) have mixed types. Specify dtype option on import or set low_memory=False.

    #beads_nhits = pd.read_csv(gzip.open(args.beads_nhits), index_col=0, squeeze=True)
    #samps_nhits = pd.read_csv(gzip.open(args.samps_nhits), index_col=0, squeeze=True)
    #hits[(beads_nhits > 2) | (samps_nhits < 2)] = 0

    columns = ['id', 'Species', 'Organism', 'Entry', 'peptide']
    hits2 = pd.merge(lib[columns], hits.reset_index(), on='id').set_index(columns)
    #print(hits.head())
    #print(hits2[hits.columns[0]])

    #print(hits2.head())
    #                                                                                                                                       S148
    #id Species               Organism                                           Entry  peptide                                                  
    #1  Papiine herpesvirus 2 Cercopithecine herpesvirus 16 (CeHV-16) (Herpes... A0A126 MRSLLFVVGAWVAALVTNLTPDAALASGTTTTAAAGNTSATASPGDN...  False
    #2  Papiine herpesvirus 2 Cercopithecine herpesvirus 16 (CeHV-16) (Herpes... A0A126 TTTTAAAGNTSATASPGDNATSIDAGSTITAAAPPGHSTPWPALPTD...  False
    #3  Papiine herpesvirus 2 Cercopithecine herpesvirus 16 (CeHV-16) (Herpes... A0A126 ITAAAPPGHSTPWPALPTDLALPLVIGGLCALTLAAMGAGALLHRCC...  False
    #4  Papiine herpesvirus 2 Cercopithecine herpesvirus 16 (CeHV-16) (Herpes... A0A126 LCALTLAAMGAGALLHRCCRRCARRRQNVSSVSA                  False
    #5  Papiine herpesvirus 2 Cercopithecine herpesvirus 16 (CeHV-16) (Herpes... A0A130 RDRGPSRSRVRYTRLAASEA                                False


    # virus_scores = calc_virus_scores(hits2[hits.name], args.level, args.epitope_len, args.max_mismatch)
    #virus_scores = calc_virus_scores(hits2[hits.name], args.level, args.epitope_len)
    virus_scores = calc_virus_scores(hits2[hits.columns[0]], args.level, args.epitope_len)
    virus_scores.to_csv(sys.stdout, header=True)


