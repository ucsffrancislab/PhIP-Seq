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
#
#	looks like minimum l is 4


#	The resulting files may contain peptides that are identical or highly similar to each other. Eliminate some of this redundancy using the cd-hit tool, similar to what is done in the UniProt database, using the following commands:
#
#	cd-hit -i orf_tiles.fasta -o orf_tiles_clustered.fasta -c 0.95 -G 0 -A 50 -M 0 -T 1 -d 0
#
#	cd-hit -i cterm_tiles.fasta -o cterm_tiles_clustered.fasta -c 0.95 -G 0 -aL 1.0 -aS 1.0 -M 0 -T 1 -d 0
#
#	In this particular case, we are clustering the peptide tiles to 95% (-c 0.95) local identity (-G 0) while controlling the alignment coverage (-A 50 requires the alignment to cover at least 50 amino acids). The C-terminal peptides are aligned more stringently to ensure that the final residues of the ORF are not lost (-aL 1.0 –aS 1.0 requires 100% of each sequence to be aligned with possible mismatches). Specifying –M 0 allows unlimited memory, -T 1 specifies one CPU thread, and –d 0 ensures sequence names are not truncated. See cd-hit documentation for more options (cd-hit.org). The clustered peptides are written to orf_tiles_clustered.fasta and cterm_tiles_clustered.fasta.

#	-c	sequence identity threshold, default 0.9
#		this is the default cd-hit's "global sequence identity" calculated as:
#		number of identical amino acids or bases in alignment
#		divided by the full length of the shorter sequence
#	-G	use global sequence identity, default 1
#		if set to 0, then use local sequence identity, calculated as :
#		number of identical amino acids or bases in alignment
#		divided by the length of the alignment
#		NOTE!!! don't use -G 0 unless you use alignment coverage controls
#		see options -aL, -AL, -aS, -AS
#	-b	band_width of alignment, default 20
#	-T	number of threads, default 1; with 0, all CPUs will be used
#	-M	memory limit (in MB) for the program, default 800; 0 for unlimitted;
#	-l	length of throw_away_sequences, default 10					# <===== minimum 4
#	-n	word_length, default 5, see user's guide for choosing it
#	-t	tolerance for redundance, default 2
#	-d	length of description in .clstr file, default 20
#		if set to 0, it takes the fasta defline and stops at first space
#	-s	length difference cutoff, default 0.0
#		if set to 0.9, the shorter sequences need to be
#		at least 90% length of the representative of the cluster
#	-S	length difference cutoff in amino acid, default 999999
#		if set to 60, the length difference between the shorter sequences
#		and the representative of the cluster can not be bigger than 60
#	-aL	alignment coverage for the longer sequence, default 0.0
#		if set to 0.9, the alignment must covers 90% of the sequence
#	-AL	alignment coverage control for the longer sequence, default 99999999
#		if set to 60, and the length of the sequence is 400,
#		then the alignment must be >= 340 (400-60) residues
#	-aS	alignment coverage for the shorter sequence, default 0.0
#		if set to 0.9, the alignment must covers 90% of the sequence
#	-AS	alignment coverage control for the shorter sequence, default 99999999
#		if set to 60, and the length of the sequence is 400,
#		then the alignment must be >= 340 (400-60) residues
#	-A	minimal alignment coverage control for the both sequences, default 0
#		alignment must cover >= this value for both sequences 
#	-uL	maximum unmatched percentage for the longer sequence, default 1.0
#		if set to 0.1, the unmatched region (excluding leading and tailing gaps)
#		must not be more than 10% of the sequence
#	-uS	maximum unmatched percentage for the shorter sequence, default 1.0
#		if set to 0.1, the unmatched region (excluding leading and tailing gaps)
#		must not be more than 10% of the sequence
#	-U	maximum unmatched length, default 99999999
#		if set to 10, the unmatched region (excluding leading and tailing gaps)
#		must not be more than 10 bases
#	-B	1 or 0, default 0, by default, sequences are stored in RAM
#		if set to 1, sequence are stored on hard drive
#		!! No longer supported !!
#	-p	1 or 0, default 0
#		if set to 1, print alignment overlap in .clstr file
#	-g	1 or 0, default 0
#		by cd-hit's default algorithm, a sequence is clustered to the first 
#		cluster that meet the threshold (fast cluster). If set to 1, the program
#		will cluster it into the most similar cluster that meet the threshold
#		(accurate but slow mode)
#		but either 1 or 0 won't change the representatives of final clusters
#	-sc	sort clusters by size (number of sequences), default 0, output clusters by decreasing length
#		if set to 1, output clusters by decreasing size
#	-sf	sort fasta/fastq by cluster size (number of sequences), default 0, no sorting
#		if set to 1, output sequences by decreasing cluster size
#		this can be very slow if the input is in .gz format


#	~/.local/cd-hit-v4.8.1-2019-0228/cd-hit \
#		-i orf_tiles-${TILESIZE}-${OVERLAP}.fasta \
#		-o orf_tiles_clustered-${TILESIZE}-${OVERLAP}.fasta \
#		-c 1.00 -G 0 -A 100 -M 0 -T 1 -d 0 -l 4
#	#	-c 0.95 -G 0 -A  50 -M 0 -T 1 -d 0 -l 4	# <-- Try this again?
#	
#	#	-c 1.00 -G 0 -A 100 -M 0 -T 1 -d 0 -l $((TILESIZE/10+4))
#	#	-c 0.95 -G 0 -A 50 -M 0 -T 1 -d 0 -l $((TILESIZE/10+4))	#	from paper
#	#	-c 0.95 -G 0 -A 100 -M 0 -T 1 -d 0 -l $((TILESIZE/10+4))
#	#	no difference with 56/28
#	#	-c 0.95 -G 0 -A 50 -M 0 -T 1 -d 0
#
#	#	>EGFR_HotSpot-EGFR:1016-1026:Mutation:Q1021*|CTERM|STOP
#	#	YLIPQ*GFFSS
#	#	Warning: from file "cterm_tiles-56-28.fasta",
#	#	Discarding invalid sequence or sequence without identifier and description!
#	
#	~/.local/cd-hit-v4.8.1-2019-0228/cd-hit \
#		-i cterm_tiles-${TILESIZE}-${OVERLAP}.fasta \
#		-o cterm_tiles_clustered-${TILESIZE}-${OVERLAP}.fasta \
#		-c 1.00 -G 0 -A 100 -aL 1.0 -aS 1.0 -M 0 -T 1 -d 0 -l 4
#	#	-c 0.95 -G 0 -aL 1.0 -aS 1.0 -M 0 -T 1 -d 0 -l 4 #	<-- Try this again?
#	
#	#	-c 1.00 -G 0 -A 100 -aL 1.0 -aS 1.0 -M 0 -T 1 -d 0 -l $((TILESIZE/10+4))
#	#	-c 0.95 -G 0 -aL 1.0 -aS 1.0 -M 0 -T 1 -d 0 -l $((TILESIZE/10+4))	#	from paper
#	
#	
#	cat orf_tiles_clustered-${TILESIZE}-${OVERLAP}.fasta \
#		cterm_tiles_clustered-${TILESIZE}-${OVERLAP}.fasta \
#		| pepsyn pad -l $TILESIZE --c-term - - \
#		> protein_tiles-${TILESIZE}-${OVERLAP}.fasta

cat orf_tiles-${TILESIZE}-${OVERLAP}.fasta \
	cterm_tiles-${TILESIZE}-${OVERLAP}.fasta \
	| pepsyn pad -l $TILESIZE --c-term - - \
	> protein_tiles-${TILESIZE}-${OVERLAP}.fasta


PREFIX=AGGAATTCCGCTGCGT
SUFFIX=GCCTGGAGACGCCATC
PREFIXLEN=${#PREFIX}
SUFFIXLEN=${#SUFFIX}
FREQTHRESH=0.01

#	pepsyn does not seem to have ever actually implemented the --codon-table option

cat protein_tiles-${TILESIZE}-${OVERLAP}.fasta \
	| pepsyn revtrans --codon-freq-threshold $FREQTHRESH --amber-only - - \
	| pepsyn prefix -p $PREFIX - - \
	| pepsyn suffix -s $SUFFIX - - \
	| pepsyn recodesite --site EcoRI --site HindIII --clip-left $PREFIXLEN \
	--clip-right $SUFFIXLEN --codon-freq-threshold $FREQTHRESH \
	--amber-only - oligos-${TILESIZE}-${OVERLAP}.fasta


#	findsite is not specific to frames so the following are found MANY times

#	must clip-left as the prefix contain EcoRI
pepsyn findsite --site EcoRI --clip-left 3 oligos-${TILESIZE}-${OVERLAP}.fasta

pepsyn findsite --site HindIII oligos-${TILESIZE}-${OVERLAP}.fasta


pepsyn clip \
	--left $PREFIXLEN \
	--right $SUFFIXLEN \
	oligos-${TILESIZE}-${OVERLAP}.fasta \
	oligos-ref-${TILESIZE}-${OVERLAP}.fasta


#	bowtie-build -q oligos-ref.fasta bowtie_index/mylibrary





