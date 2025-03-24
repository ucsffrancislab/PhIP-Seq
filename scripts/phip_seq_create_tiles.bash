#!/usr/bin/env bash

set -e  # exit if any command fails
set -u  # Error on usage of unset variables
set -o pipefail
if [ -n "$( declare -F module )" ] ; then
  echo "Loading required modules"
#  module load CBI samtools r
fi
#set -x # print expanded command before executing it
#set -v


function usage(){
	set +x
cat <<EOF
Usage:

$0 -t 56 -o 28 -i Human_alphaherpesvirus_3_proteins.faa

	tilesize = integer

	overlap = integer less than tilesize

	input = an amino acid fasta


	PhIP-Seq Characterization of Serum Antibodies Using Oligonucleotide Encoded Peptidomes

	https://pmc.ncbi.nlm.nih.gov/articles/PMC6568263

	python3 -m pip uninstall -y biopython pepsyn
	python3 -m pip install --upgrade --user numpy scipy biopython click tqdm phip-stat
	git clone https://github.com/lasersonlab/pepsyn.git
	cd pepsyn
	python3 setup.py install --user
	cd ..

	http://weizhong-cluster.ucsd.edu/cd-hit/
	https://github.com/weizhongli/cdhit/releases
	wget https://github.com/weizhongli/cdhit/releases/download/V4.8.1/cd-hit-v4.8.1-2019-0228.tar.gz

EOF

	exit
}

while [ $# -gt 0 ] ; do
	case $1 in
		-t|--tilesize)
			shift; TILESIZE=$1; shift;;
		-o|--overlap)
			shift; OVERLAP=$1; shift;;
		-i|--input)
			shift; INPUT=$1; shift;;
		*)
			echo "Unknown param :${1}:"; usage ;;
	esac
done

#	TILESIZE=56
#	OVERLAP=28

#	TILESIZE=$1
#	OVERLAP=$2
#	INPUT=$3			#	is this a fasta file? list of peptides? what? looks like its an amino acid fasta


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




#	-l default is 10 so any read shorter is ignored
#		set to ..... half of tile size????
#	$((TILESIZE/10))

#	looks like minimum l is 4


~/.local/cd-hit-v4.8.1-2019-0228/cd-hit \
	-i orf_tiles-${TILESIZE}-${OVERLAP}.fasta \
	-o orf_tiles_clustered-${TILESIZE}-${OVERLAP}.fasta \
	-c 0.95 -G 0 -A 100 -M 0 -T 1 -d 0 -l $((TILESIZE/10+4))
#	no difference with 56/28
#	-c 0.95 -G 0 -A 50 -M 0 -T 1 -d 0

~/.local/cd-hit-v4.8.1-2019-0228/cd-hit \
	-i cterm_tiles-${TILESIZE}-${OVERLAP}.fasta \
	-o cterm_tiles_clustered-${TILESIZE}-${OVERLAP}.fasta \
	-c 0.95 -G 0 -aL 1.0 -aS 1.0 -M 0 -T 1 -d 0 -l $((TILESIZE/10+4))


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





