#!/usr/bin/env python3
from __future__ import division

import os
import sys
#import gzip
import difflib
import argparse

# import regex
import pandas as pd



#	This is WAY, WAY FASTER. Got it working just as I didn't need it anymore.
class PeptideIndex:
	def __init__(self, epitope_len):
		self.epitope_len = epitope_len
		self.kmers = set()

	def add_peptide(self, peptide):
		# Extract all substrings of length epitope_len
		for i in range(len(peptide) - self.epitope_len + 1):
			self.kmers.add(peptide[i:i + self.epitope_len])

	def has_match(self, peptide):
		# Check if any substring exists in our index
		for i in range(len(peptide) - self.epitope_len + 1):
			if peptide[i:i + self.epitope_len] in self.kmers:
				return True
		return False


##	This is 3-6 times faster
#def is_novel_peptide_fast(peptide, assigned_peptides, epitope_len):
#	# Get all kmers from the query peptide
#	query_kmers = {peptide[i:i + epitope_len]
#		for i in range(len(peptide) - epitope_len + 1)}
#
#	# Check against all assigned peptides
#	for assigned in assigned_peptides:
#		for i in range(len(assigned) - epitope_len + 1):
#			if assigned[i:i + epitope_len] in query_kmers:
#		return False
#	return True
#
##
##	This will get progressively slower and slower and the list of assigned peptides gets longer and longer
##
#def is_novel_peptide(peptide, assigned_peptides, epitope_len):
#	for assigned_peptide in assigned_peptides:
#		matcher = difflib.SequenceMatcher(None, peptide, assigned_peptide)
#		match = matcher.find_longest_match(0, len(peptide), 0, len(assigned_peptide))
#		#	the actual length is irrelevant
#		if match.size >= epitope_len:
#			return False
#	return True
#
#def is_novel_peptide(peptide, assigned_peptides, epitope_len, max_mm):
#	for i in range(len(peptide) - epitope_len):
#		epitope = peptide[i:i + epitope_len]
#		fmt = r'({epitope}){{s<={max_mm}}}'.format(epitope=epitope, max_mm=max_mm)
#		r = regex.compile(fmt)
#		for assigned_peptide in assigned_peptides:
#			if r.search(assigned_peptide) is not None:
#				return False
#		return True


def calc_virus_scores(series, level, epitope_len, species_order):
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

	index = PeptideIndex(epitope_len)

	assigned_peptides = set()
	peptides_for_export = []
	grouped = series[series].groupby(level=level)
	#for virus in nhits_per_virus[nhits_per_virus > 0].order(ascending=False).index:

	#print(nhits_per_virus[nhits_per_virus > 0].sort_values(ascending=False))
	#print(nhits_per_virus[nhits_per_virus > 0].sort_values(ascending=False).head())
	#Species
	#Influenza A virus           100
	#Human herpesvirus 4          84
	#Enterovirus B                44
	#Rhinovirus A                 24
	#Streptococcus pneumoniae     23
	#Name: S148, dtype: int64

#    for virus in nhits_per_virus[nhits_per_virus > 0].sort_values(ascending=False).index:

	for virus in species_order:

		if virus not in grouped.groups:
			continue

		virus_hits = grouped.get_group(virus)

		#   These the True tiles

		#print(type(virus_hits))
		#<class 'pandas.core.series.Series'>

		#print(virus_hits.index[0:4])
		#MultiIndex([(  42, 'Human herpesvirus 5', ...),
		#    (2858, 'Human herpesvirus 5', ...),
		#    (2870, 'Human herpesvirus 5', ...),
		#    (2905, 'Human herpesvirus 5', ...)],
		#   names=['id', 'Species', 'peptide'])

		#print(virus_hits.head())
		# id    Species            peptide
		# 359   Influenza A virus  WYGYHHQNEQGSGYAADKKSTQNAIDGITNKVNSVIEKMNTQFTAVGKEFNKLERR    True
		# 535   Influenza A virus  HGLKRGPSTEGVPESMREEYRKEQQSAVDADDSHFVNIELE                   True
		# 553   Influenza A virus  IGNGCFEFYHKCDNECMESVRNGTYDYPKYSEESKLNREKIDGVKLESMGVYQILA    True
		# 1133  Influenza A virus  DAPFLDRLRRDQKSLKGRGSTLGLNIETATLAGKQIVEWILKEESNETLKMSIASV    True
		# 1135  Influenza A virus  PSSRYLADMTLEEMSRDWFMLMPRQKVAGSLCVRMDQAIMEKNIVLKANFSVIFDR    True
		# Name: S1, dtype: bool

		score = 0
		#peptides = virus_hits.index.get_level_values('peptide')

		#print(peptides.shape)
		# (341,)

		#print(peptides)
		# Index(['WYGYHHQNEQGSGYAADKKSTQNAIDGITNKVNSVIEKMNTQFTAVGKEFNKLERR',
		#        'HGLKRGPSTEGVPESMREEYRKEQQSAVDADDSHFVNIELE',
		#        'IGNGCFEFYHKCDNECMESVRNGTYDYPKYSEESKLNREKIDGVKLESMGVYQILA',
		#        'DAPFLDRLRRDQKSLKGRGSTLGLNIETATLAGKQIVEWILKEESNETLKMSIASV',
		#        'PSSRYLADMTLEEMSRDWFMLMPRQKVAGSLCVRMDQAIMEKNIVLKANFSVIFDR',
		#        'IDLADSEMNKLYERVKRQLRENAEEDGTGCFEIFHKCDDDCMASIRNNTYDHSKYR',
		#        'DIWVTREPYVSCDTSKCYQFALGQGTTLNNKHSN',
		#        'DGWYGYHHSNEQGSGYAADQESTQKAIDGVTNKVNSIINKMNTQFEAVGREFNNLE',


		# I feel like the assigned_peptides needs to be reset otherwise it compounds.
		# This would change the resulting virus score based on the order in which it were processed.
		# The viruses are processed sorted from highest to lowest count.
		#assigned_peptides.clear()

		#print(len(assigned_peptides))

		#print(peptides[0:2])
		#Index(['PGDTIIFEANGNLIAPWYAFALSRGFGSGIITSNASMGECDAKCQTPQGAINSSLP',
		#  'QIASNENMETIDSITLELRSKYWAIRTRSGGNTNKQRASAGQISVQPTFSVQRNLP'], dtype='object', name='peptide')

		#for peptide in peptides:
		for virus_hit in virus_hits.index:
			#print(virus_hit)
			#(42, 'Human herpesvirus 5', 'TIKNTKPQCRPEDYATRLQDLRVTFHRVKPTLQREDDYSVWLDGDHVYPGLKTELH')
			peptide=virus_hit[2]


			#if is_novel_peptide(peptide, assigned_peptides, epitope_len, max_mm):
			#if is_novel_peptide(peptide, assigned_peptides, epitope_len):
			#if is_novel_peptide_fast(peptide, assigned_peptides, epitope_len):
			if not index.has_match(peptide):
				score += 1
				assigned_peptides.add(peptide)
				index.add_peptide(peptide)
				peptides_for_export.append([virus_hit[0],virus_hit[1],virus_hit[2]])

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

		#print(virus + " : " + str(len(assigned_peptides)))
		virus_scores[virus] = score

	pd.DataFrame(peptides_for_export, columns=['id', 'species', 'peptide']).sort_values(by=['id']).to_csv(
		args.hits.rsplit('.',1)[0]+'.'+str(args.epitope_len[0])+'.peptides.txt',index=False)

	return virus_scores


if __name__ == '__main__':
#    parser = argparse.ArgumentParser()
#    parser.add_argument("hits")
#    parser.add_argument("oligo_metadata")
#    #parser.add_argument("beads_nhits")
#    #parser.add_argument("samps_nhits")
#    parser.add_argument("level")
#    parser.add_argument("epitope_len", type=int)
#    # parser.add_argument("max_mismatch", type=int)
#    args = parser.parse_args()



	parser = argparse.ArgumentParser(prog=os.path.basename(__file__))

#	parser.add_argument('files', nargs='*', help='files help')
#	parser.add_argument('-V','--version', action='version', version='%(prog)s 1.1')
#	parser.add_argument('-o', '--output', nargs=1, type=str, default=['merged.csv.gz'], help='output csv filename to %(prog)s (default: %(default)s)')
#	#parser.add_argument('-s', '--sep', nargs=1, type=str, default='\t', help='the separator to %(prog)s (default: %(default)s)')
#	#	store_true means "int=False unless --int passed, then int=True" (store_false is the inverse)
#	parser.add_argument('--int', action='store_true', help='convert values to ints to %(prog)s (default: %(default)s)')
#	parser.add_argument('--seqint', action='store_true', help='items are ints so sort like it : %(prog)s (default: %(default)s)')


# elledge_calc_scores_nofilter.py All.count.Zscores.${sample}.csv /francislab/data1/refs/PhIP-Seq/VIR3_clean.virus_score.csv Species ${length} > tmp

	parser.add_argument("--hits", required=True)
	parser.add_argument("--oligo_metadata", required=True)
	parser.add_argument("--level", default=['species'])
	parser.add_argument("--epitope_len", type=int, default=[7])
	parser.add_argument("--species_order", required=True)

	# read arguments from the command line
	args = parser.parse_args()





	species_order = data=[line.strip() for line in open(args.species_order,'r')]
	#print(species_order[0:5])
	#['Human herpesvirus 5', 'Human herpesvirus 2', 'Influenza A virus', 'Human herpesvirus 4', 'Enterovirus C']

	hits = pd.read_csv(args.hits, index_col=0)
	#print(hits.head())
	#       S1
	#id
	#1   False
	#2   False
	#3   False
	#5    True
	#6   False

	lib = pd.read_csv(args.oligo_metadata,low_memory=False)
	#print(lib.columns)
	#Index(['Unnamed: 0', 'Aclstr50', 'Bclstr50', 'Entry', 'Gene names',
	#       'Gene ontology (GO)', 'Gene ontology IDs', 'Genus', 'Organism',
	#       'Protein names', 'Sequence', 'Species', 'Subcellular location',
	#       'Version (entry)', 'Version (sequence)', 'end', 'id', 'oligo', 'source',
	#       'start', 'peptide'],
	#      dtype='object')
	# DtypeWarning: Columns (7) have mixed types. Specify dtype option on import or set low_memory=False.

	#beads_nhits = pd.read_csv(gzip.open(args.beads_nhits), index_col=0, squeeze=True)
	#samps_nhits = pd.read_csv(gzip.open(args.samps_nhits), index_col=0, squeeze=True)
	#hits[(beads_nhits > 2) | (samps_nhits < 2)] = 0

	#columns = ['id', 'Species', 'Organism', 'Entry', 'peptide']
	#columns = ['id', 'Species', 'peptide']
	columns = ['id', 'species', 'peptide']
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

	virus_scores = calc_virus_scores(hits2[hits.columns[0]], args.level[0], args.epitope_len[0], species_order)
	virus_scores.to_csv(sys.stdout, header=True)



