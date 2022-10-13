#!/bin/zsh
#SBATCH --job-name=FillTaxids
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2g
#SBATCH --time=100:00:00
#SBATCH --account=tyjames1
#SBATCH --partition=standard

#Blast output filei (IN) and path to blast database (DB)
IN=$1
DB=$2

#Check arguments
if [[ -z $IN || -z $DB ]]; then
    echo "ERROR: Expected two positional arguments."
    echo "USAGE: sh 2_uniref_blastout_get_and_fill.sh [blast_output] [database]"
    exit 0
fi
    
#blastdbcmd output file to write
TEMP=${IN}.blastdbcmd.out.tmp

#Output file with subject accessions and taxids
TAXIDS=${IN}.mappable_taxids

#Final blast output with taxids filled-in
OUT=${IN}.TaxIDified

#Database that was blasted against
DATABASE=${DB}

IFS='\t'
while read p; do
    subject_accession=`echo $p | cut -f2`
    echo $subject_accession
    blastdbcmd -entry "$subject_accession" -db "$DATABASE" -outfmt "%a;%t" >> $TEMP
done < $IN
unset IFS

cat $TEMP | sed -r "s/^(.*);.*TaxID=([a-zA-Z0-9_]*).*$/\1;\2/" > $TAXIDS

#rm $TEMP

python ~/scripts/alien_index_scores/util/blastout_fill_taxids.py -b $IN -t $TAXIDS > ${OUT} 
