#!/usr/bin/env bash

# Directory to work in. Put in the absolute path (starting with /).
WORKDIR="$HOME/test/BLAST"

# Name of the fasta file in the working directory that contains the query sequence.
QUERY_SEQ="Bacillus_Cereus.fsa"

# Name of the csv file in the working directory that has the GenBank accession
# numbers in the third column. To change the behavior, you'll have to modify
# line 69 (the line with the awk command).
ACCESSION_NUMS="Bacillus_strains.csv"

# A human-readable name for the constructed BLAST database.
BLAST_DB_NAME="22_Bacillus_strains"

# Whether or not to verify the locally built BLAST database (could be very long output).
VERIFY_BLAST_DB=true


##############################################################################
# Shouldn't need to modify below this line
##############################################################################

function is_exe() {
  which "$@" &> /dev/null
}

function setup() {
  ## Create directories for analysis
  echo "Changing directory to ${WORKDIR}"
  mkdir -p "${WORKDIR}" && cd "${WORKDIR}" || exit 1
  mkdir -p blastdb queries fasta results blastdb_custom

  ## Add local software path to $PATH
  if ! is_exe conda; then
    echo "It seems Anaconda isn't installed." \
    "Follow instructions on https://docs.conda.io/en/latest/miniconda.html" \
    "to install it to your system."
  fi

  ## Check the executables used
  for cmd in esearch elink efetch makeblastdb blastdbcmd blastn; do
    if ! is_exe $cmd; then
      echo "Error: ${cmd} not found in ${PATH}"
      echo "Try this command: conda install -c bioconda entrez-direct blast"
      echo "If it's already installed, maybe you didn't activate the conda environment?"
      exit 1
    fi
  done
}

function get_query_seq() {
  if [ -f "$QUERY_SEQ" ]; then
    echo "  Copying query sequence ${QUERY_SEQ}"
    cp "${WORKDIR}/${QUERY_SEQ}" "queries/${QUERY_SEQ}"
  else
    echo "Error: Query sequence ${QUERY_SEQ} not found in ${WORKDIR}"
    exit 1
  fi
}

function get_db_seq() {
  if [ ! -f "$ACCESSION_NUMS" ]; then
    echo "Error: csv file for accession numbers (${ACCESSION_NUMS}) not found in ${WORKDIR}"
    exit 1
  fi
  ## Get GenBank accession numbers
  strains=$(awk -F, 'BEGIN {ORS=" OR "} NR>1 {print $3}' "${ACCESSION_NUMS}" | paste -d, -s | sed 's/ OR $//')

  ## Skip download when fasta file is found in the target directory.
  ## Otherwise download files using the Entrez Direct suite.
  if [ -f "fasta/${BLAST_DB_NAME}.fsa" ]; then
    echo "Sequences found in ${WORKDIR}/fasta/${BLAST_DB_NAME}.fsa, skipping download."
  else
    echo "  Retriving database sequences using ${ACCESSION_NUMS}"
    echo "Query: ${strains}"
    echo "    This could take a while..."

    esearch \
        -db assembly \
        -query "${strains}" |
    elink \
        -target nucleotide \
        -name assembly_nuccore_insdc |
    efetch \
        -format fasta > "fasta/${BLAST_DB_NAME}.fsa"

    printf "\n\nSaved downloaded sequences in fasta/%s.fsa\n" "${BLAST_DB_NAME}"
  fi
}

function gen_taxid_map() {
  if [ -f taxid_map.txt ]; then
    echo "taxid_map file found"
  else
    echo "Preparing taxid_map file"
    touch taxid_map.txt
    for i in $(echo $strains | sed 's/ OR / /g'); do
      echo $i
      line=$(esearch -db assembly -query $i | elink -target taxonomy | efetch -format uid)
      echo "${i} ${line}"  >> taxid_map.txt
    done
  fi
}

function make_blast_db() {
  makeblastdb \
      -in "${WORKDIR}/fasta/${BLAST_DB_NAME}.fsa" \
      -dbtype nucl \
      -parse_seqids -out $BLAST_DB_NAME \
      -title $BLAST_DB_NAME \
      -taxid_map "${WORKDIR}/taxid_map.txt" \
      -blastdb_version 5
}

function verify_blast_db() {
  printf "\n\n Accession     Sequence length\n"
  blastdbcmd \
      -entry all \
      -db $BLAST_DB_NAME \
      -outfmt "%a %l"
}

function run_blastn() {
  blastn \
      -query "${WORKDIR}/queries/${QUERY_SEQ}" \
      -db "${WORKDIR}/blastdb_custom/${BLAST_DB_NAME}" \
      -html \
      -out "${WORKDIR}/results/${QUERY_SEQ}.html"
  blastn \
      -query "${WORKDIR}/queries/${QUERY_SEQ}" \
      -db "${WORKDIR}/blastdb_custom/${BLAST_DB_NAME}" \
      -out "${WORKDIR}/results/${QUERY_SEQ}.out"
}

##############################################################################
# Step 0. Setup
##############################################################################
## Create working directories, and check for software used
setup

##############################################################################
# Step 1. Retrieve sequences
##############################################################################
echo "Step 1: Retrieve sequences"

## Retrieve query sequence and copy it to the right directory
get_query_seq

## Retrieve database sequences
get_db_seq

##############################################################################
# Step 2. Make BLAST database
##############################################################################
echo "Step 2: make BLAST database"

## Generate taxid_map file
gen_taxid_map

## Make local database
cd $WORKDIR/blastdb_custom || exit
make_blast_db

## Verify BLAST database
if [ $VERIFY_BLAST_DB = true ]; then
  verify_blast_db
fi

## Go back to working directory
cd $WORKDIR || exit

##############################################################################
# Step 3. Run BLAST+
##############################################################################
echo "Step 3: run BLAST+"
run_blastn && \
  echo "Saved results to ${WORKDIR}/results/${QUERY_SEQ}.html"
