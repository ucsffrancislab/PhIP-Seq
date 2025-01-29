#!/usr/bin/env bash



#	pepsyn is a bit outdated. Need biopython with Alphabet (1.69)
#	but I can't get it to install anymore.
#	I think that there is a different manner of creating these.

#	https://pmc.ncbi.nlm.nih.gov/articles/PMC6568263


#	pip uninstall biopython pepsyn
#	python3 -m pip install --upgrade --user numpy scipy biopython click tqdm phip-stat
#	cd ~/github/
#	mkdir lasersonlab
#	cd lasersonlab/
#	git clone https://github.com/lasersonlab/pepsyn.git
#	cd pepsyn
#	python3 setup.py install --user



#	python3 -m pip install --upgrade --user numpy scipy biopython==1.69 click tqdm pepsyn phip-stat
#	
#	http://weizhong-cluster.ucsd.edu/cd-hit/
#	https://github.com/weizhongli/cdhit/releases
#	wget https://github.com/weizhongli/cdhit/releases/download/V4.8.1/cd-hit-v4.8.1-2019-0228.tar.gz


#	TILESIZE=56
#	OVERLAP=28

TILESIZE=$1
OVERLAP=$2
INPUT=$3			#	is this a fasta file? list of peptides? what? looks like its a fasta

#	./phipseq.bash 56 28 Human_alphaherpesvirus_3_proteins.fa > phipseq.log

#head /francislab/data1/working/20230413-VZV/Human_alphaherpesvirus_3_proteins.fa 
#>NP_040124.1_membrane_protein_V1_Human_alphaherpesvirus_3
#MSRVSEYGVPEGVRESDSDTDSVFMYQHTELMQNNASPLVVQTRPPAVLI
#PLVDVPRPRSRRKASAQLKMQMDRLCNVLGVVLQMATLALVTYIAFVVHT
#RATSCKRE
#>NP_040125.2_myristylated_tegument_protein_CIRC_Human_alphaherpesvirus_3
#MGSTLVRPSLNATAEENPASETRCLLRVLAGRTVDLPGGGTLHITCTKTY
#VIIGKYSKPGERLSLARLIGRAMTPGGARTFIILAMKEKRSTTLGYECGT
#GLHLLAPSMGTFLRTHGLSNRDLCLWRGNIYDMHMQRLMFWENIAQNTTE
#TPCITSTLTCNLTEDSGEAALTTSDRPTLPTLTAQGRPTVSNIRGILKGS
#PRQQPVCHRVRFAEPTEGVLM




cat ${INPUT} \
	| pepsyn x2ggsg - - \
	| pepsyn tile -l $TILESIZE -p $OVERLAP - - \
	| pepsyn disambiguateaa - - \
	> orf_tiles-${TILESIZE}-${OVERLAP}.fasta

cat ${INPUT} \
	| pepsyn x2ggsg - - \
	| pepsyn ctermpep -l $TILESIZE --add-stop - - \
	| pepsyn disambiguateaa - - \
	> cterm_tiles-${TILESIZE}-${OVERLAP}.fasta

~/.local/cd-hit-v4.8.1-2019-0228/cd-hit \
	-i orf_tiles-${TILESIZE}-${OVERLAP}.fasta \
	-o orf_tiles_clustered-${TILESIZE}-${OVERLAP}.fasta \
	-c 0.95 -G 0 -A 100 -M 0 -T 1 -d 0
#	no difference with 56/28
#	-c 0.95 -G 0 -A 50 -M 0 -T 1 -d 0

~/.local/cd-hit-v4.8.1-2019-0228/cd-hit \
	-i cterm_tiles-${TILESIZE}-${OVERLAP}.fasta \
	-o cterm_tiles_clustered-${TILESIZE}-${OVERLAP}.fasta \
	-c 0.95 -G 0 -aL 1.0 -aS 1.0 -M 0 -T 1 -d 0


cat orf_tiles_clustered-${TILESIZE}-${OVERLAP}.fasta \
	cterm_tiles_clustered-${TILESIZE}-${OVERLAP}.fasta \
	| pepsyn pad -l $TILESIZE --c-term - - \
	> protein_tiles-${TILESIZE}-${OVERLAP}.fasta


PREFIX=AGGAATTCCGCTGCGT
SUFFIX=GCCTGGAGACGCCATC
PREFIXLEN=${#PREFIX}
SUFFIXLEN=${#SUFFIX}
FREQTHRESH=0.01

cat protein_tiles-${TILESIZE}-${OVERLAP}.fasta \
	| pepsyn revtrans --codon-freq-threshold $FREQTHRESH --amber-only - - \
	| pepsyn prefix -p $PREFIX - - \
	| pepsyn suffix -s $SUFFIX - - \
	| pepsyn recodesite --site EcoRI --site HindIII --clip-left $PREFIXLEN \
	--clip-right $SUFFIXLEN --codon-freq-threshold $FREQTHRESH \
	--amber-only - - \
	> oligos-${TILESIZE}-${OVERLAP}.fasta


pepsyn findsite --site EcoRI --clip-left 3 oligos-${TILESIZE}-${OVERLAP}.fasta

pepsyn findsite --site HindIII oligos-${TILESIZE}-${OVERLAP}.fasta


pepsyn clip \
	--left $PREFIXLEN \
	--right $SUFFIXLEN \
	oligos-${TILESIZE}-${OVERLAP}.fasta \
	oligos-ref-${TILESIZE}-${OVERLAP}.fasta


#	bowtie-build -q oligos-ref.fasta bowtie_index/mylibrary





