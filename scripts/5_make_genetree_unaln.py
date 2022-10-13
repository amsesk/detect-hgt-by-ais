import sys
import argparse
import os
import subprocess
import shutil
from scgid.sequence import DNASequenceCollection
from scgid.parsers import BlastoutParser
import logging
EVALUE_CUTOFF="1e-50"
OUTFMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids"

def index_fasta(fasta_path):
    cmd = ['samtools',
            'faidx',
            fasta_path]
    p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    out,err = p.communicate()
    if p.returncode == 0:
        return None
    else:
        logging.critical(f"samtools indexing failed with error: {err}")
        sys.exit(1)

def make_directory(path):
    if os.path.isdir(path):
        shutil.rmtree(path)

    os.mkdir(path)

    return os.path.abspath(path)

def segregate_taxrpt(target, taxrpt_path) -> list:
    pertinent_lines = []
    with open(taxrpt_path) as taxrpt:
        for line in taxrpt:
            spl = [x.strip() for x in line.split('\t')]
            query = spl[0]
            if query == target:
                pertinent_lines.append(spl)
    return pertinent_lines

def write_tsv_from_list(tsv_list, outpath):
    with open(outpath, 'w') as out:
        for line in tsv_list:
            if isinstance(line, list):
                out.write('\t'.join(line))
                out.write('\n')
            elif isinstance(line, str):
                out.write(line)
                out.write('\n')
            else:
                raise TypeError
    return None

def pull_subject_contigs_from_database(wanted_path, fasta_path, fasta_out_path,  use_grabrust) -> DNASequenceCollection:
    if use_grabrust:
        cmd = ['grabrust_these_contigs',
                '--fasta', fasta_path,
                '--wanted', wanted_path]

        p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        out,err = p.communicate()

        if p.returncode == 0:
            with open(fasta_out_path, 'w') as fasta_out:
                fasta_out.write(out.decode('utf-8'))
            subject_fasta = DNASequenceCollection().from_fasta(fasta_out_path)
            return subject_fasta

        else:
            raise RuntimeError(err)

    else:
        # Use grab_these_contigs.py instead, UNIMPLEMENTED
        pass

def blast_query_against_other_genomes(query_fasta_path, others_tsv, outdir):
    make_directory(os.path.join(outdir, "other_blastouts"))
    with open(others_tsv, 'r') as otsv:
        blastouts = []
        for line in otsv:
            spl = [x.strip() for x in line.split('\t')]
            strain = spl[0]
            proteome_path = spl[1]

            ### Make a blast database for the other proteome if it doesn't exist
            if not os.path.isfile(f"{proteome_path}.phr"):
                make_blastdb(
                        fasta = proteome_path,
                        title = os.path.basename(proteome_path),
                        outpath = proteome_path)

            ### BLAST it
            blastout_outpath = os.path.join(outdir, "other_blastouts", f"{strain}.blastout")
            cmd = ['blastp',
                    '-query', query_fasta_path,
                    '-db', proteome_path,
                    '-evalue', EVALUE_CUTOFF,
                    '-outfmt', f"{OUTFMT}",
                    '-num_threads', '1',
                    '-max_target_seqs', '1',
                    '-out', blastout_outpath]
            p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            out, err = p.communicate()

            if p.returncode == 0:
                blastouts.append( (proteome_path, blastout_outpath) )

            else:
                logging.critical(f"BLASTing against {proteome_path} failed with error: {err}")
                sys.exit(1)
    return blastouts

def make_blastdb(fasta, title, outpath):
    cmd = ['makeblastdb',
            '-in', fasta,
            '-dbtype', 'prot',
            '-title', title,
            '-out', outpath]
    p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    out,err = p.communicate()
    if p.returncode == 0:
        logging.info(f"Made BLAST database for {title}.")
        return None

    else:
        logging.critical(f"makeblastdb for {title} failed with error: {err}")
        sys.exit(1)

def blastouts_to_best_fastas (blastouts):
    other_sequences_to_add = []
    for proteome_path, blastout_path in blastouts:

        if not os.path.isfile(f"{proteome_path}.fai"):
            index_fasta(proteome_path)

        bop = BlastoutParser()
        bop.load_from_file(blastout_path)
        bop.get_best_hits()

        assert len(bop.best_hits.values()) <= 1, "More than one best blast hit. This shouldn't happen"

        best_hits = list(bop.best_hits.values())

        if len(best_hits) != 0:
            ### Write best hits
            write_tsv_from_list(bop.best_hits, f"{blastout_path}.best")

            ### Write best hit header to file for grabrust
            wanted_path = f"{blastout_path}.best.header"
            write_tsv_from_list([best_hits[0][1]], wanted_path)

            ### Pull Best Subject Sequences from proteome and add to list of DNASequence objects
            best_hit_fasta_outpath =  f"{blastout_path}.besthit.fasta"
            pull_subject_contigs_from_database(
                    wanted_path = wanted_path,
                    fasta_path = proteome_path,
                    fasta_out_path = best_hit_fasta_outpath,
                    use_grabrust = True)

            besthit = list(DNASequenceCollection().from_fasta(best_hit_fasta_outpath).seqs())[0]
            other_sequences_to_add.append(besthit)

    return other_sequences_to_add

logging.basicConfig(level=logging.DEBUG)

### Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("-c", "--candidates", required = True, action = "store", help = "Path to summary of high AIS proteins generated by 4_")
parser.add_argument("-t", "--taxrpt", required = True, action = "store", help = "Path to taxrpt output generated by 3_")
parser.add_argument("-q", "--query", required = True, action = "store", help = "Path to protein FASTA for genome you're trying to detect HGT in.")
parser.add_argument("-d", "--database", required = True, action = "store", help = "Path to the database that you originally BLASTed against.")
parser.add_argument("--blast_others", required = False, action = "store", help = "Path to a two-column tsv of other genomes against which to BLAST each candidate. Column 1 should be a strain/species designation and Column 2 should be the path to each predicted proteome. Use this when you want to include other tips in the tree that aren't in the database, perhaps because you removed them to get accurate AISs.")
parser.add_argument("--grabrust", required = False, action = "store_true", help = "Use grabrust-these-contigs to get contigs from the database FASTA instead of a python implementation. This is much faster, but you need to have it built and in PATH.")
args = parser.parse_args()

candidates = args.candidates
taxrpt = args.taxrpt
query_genome = args.query
db = args.database
others_tsv = args.blast_others
use_grabrust = args.grabrust

if not os.path.isfile(f"{query_genome}.fai"):
    index_fasta(query_genome)

with open(candidates, 'r') as cand:
    for line in cand:
        spl = [x.strip() for x in line.split("\t")]
        this_candidate = spl[0]

        logging.info(f"Started working on {this_candidate}")

        cwd = make_directory(this_candidate)

        ### Generate and write segregate taxrpt
        seg_taxrpt_outpath = os.path.join(cwd, f"{this_candidate}.segregate_taxrpt")
        seg_taxrpt = segregate_taxrpt(target = this_candidate, taxrpt_path = taxrpt)
        write_tsv_from_list(seg_taxrpt, seg_taxrpt_outpath)

        ### Generate and write wanted list for grabbing
        wanted_outpath = os.path.join(cwd, f"{this_candidate}.wanted")
        subjects_to_pull = [x[1] for x in seg_taxrpt]
        logging.info(f"Found {len(subjects_to_pull)} subject sequences that hit to the {this_candidate} protein in its taxrpt.")
        write_tsv_from_list(subjects_to_pull, wanted_outpath)

        ### Pull wanted subject contigs from database fasta
        logging.info("Pulling sequences from database FASTA. This can take awhile.")
        wanted_fasta_outpath = os.path.join(cwd, f"{this_candidate}.wanted.fasta")
        wanted_fasta = pull_subject_contigs_from_database(
                wanted_path = wanted_outpath,
                fasta_path = db,
                fasta_out_path = wanted_fasta_outpath,
                use_grabrust = use_grabrust)

        logging.info(f"Recovered {len(wanted_fasta.seqs())}/{len(subjects_to_pull)} subject sequences from the database at: {db}")

        ### Pull query sequence from query genome file
        query_wanted_path = os.path.join(cwd, f'{this_candidate}.query')
        query_fasta_outpath = os.path.join(cwd, f'{this_candidate}.query.fasta')
        write_tsv_from_list([this_candidate], query_wanted_path)
        query_fasta = pull_subject_contigs_from_database(
                wanted_path = query_wanted_path,
                fasta_path = query_genome,
                fasta_out_path = query_fasta_outpath,
                use_grabrust = use_grabrust)

        ### Make sure there is only one sequence in query_fasta
        query_sequence = list(query_fasta.seqs())
        if len(query_sequence) == 0:
            logging.critical(f"Query sequence {this_candidate} not found in FASTA at: {query_genome}")
            sys.exit(1)
        elif len(query_sequence) > 1:
                logging.critical(f"More than one sequence matching query sequence in {this_candidate} in FASTA at: {query_genome}")
                sys.exit(1)
        ### Add query sequence at the head of wanted_fasta
        else:
            logging.info(f"Adding query sequence to pulled FASTA: {this_candidate}")
            unaln_sequences = list(wanted_fasta.seqs())
            unaln_sequences.insert(0, query_sequence[0])

        ### This is where we can BLAST other genomes with the query sequence, in order to add them to the unaln fasta before outputing
            if others_tsv is not None:
                blastouts = blast_query_against_other_genomes(
                        query_fasta_path = query_fasta_outpath,
                        others_tsv = others_tsv,
                        outdir = cwd)

                other_sequences_to_add = blastouts_to_best_fastas(blastouts)

                logging.info(f"Identified {len(other_sequences_to_add)} hits to {this_candidate} in other proteomes specified in {others_tsv}. Adding to unaligned fasta.")

                ### Add other sequences to unaln
                for s in other_sequences_to_add:
                    unaln_sequences.insert(1, s)

            ### Print unaligned FASTA
            unaln_outpath = os.path.join(cwd, f"{this_candidate}.unaln.fasta")
            with open(unaln_outpath, 'w') as unaln:
                for s in unaln_sequences:
                    unaln.write(f">{s.header}\n{s.string}\n")

            logging.info(f"Wrote {len(unaln_sequences)} sequences to unaligned FASTA for {this_candidate} to {unaln_outpath}")
            logging.info(f"Done with {this_candidate}.")
            logging.info("-"*25)
