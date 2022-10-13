##########################
# Written by Kevin Amses
# amsesk@umich.edu
##########################

import sys
import os
from ete3 import NCBITaxa
from collections import OrderedDict

def fill_lineage_from_supplement(taxid, supplement):
    supp_d = {}
    for line in open(supplement, 'r'):
        spl = [x.strip() for x in line.split('\t')]
        supp_d[spl[0]] = spl[1:]

    if taxid in supp_d:
        return supp_d[taxid]

    else:
        return None

if sys.argv[1] in ['-h','--help']:
    print ('''
ncbi_taxrpt.py
USAGE: python ncbi_taxrpt_all.py [blast_output] ([supplemental_lineages])\n
    REQUIREMENT IN ORIGINAL BLAST CALL:
    \t-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids"
    ''')
    sys.exit(0)
ncbi = NCBITaxa()
ids = OrderedDict()
rank_of_int = ('superkingdom','kingdom','phylum','class','order','family','genus','species')
#print '\t'.join(['query','hit','evalue'])+'\t'+'\t'.join(rank_of_int)

print ("{}\t{}".format('\t'.join(['query','hit','evalue','bitscore']), '\t'.join(rank_of_int)))

lines = []
with open(sys.argv[1], 'r') as blastout:
    for line in blastout:

        from_supplement = False

        spl = [x.strip() for x in line.split("\t")]
        name = spl[0]
        hit = spl[1]
        tid = spl[12]
        bitscore = spl[11]
        evalue = spl[10]

        if tid == "N/A":
            lineage = "No_associated_taxid"
        else:
            try:
                if len(tid.split(';')) > 1:
                    for i in tid.split(';'):
                        i = int(i)
                        lineage = ncbi.get_taxid_translator(ncbi.get_lineage(i))
                else:
                    tid = int(tid)
                    lineage = ncbi.get_taxid_translator(ncbi.get_lineage(tid))
            except:
                try:
                    lineage = fill_lineage_from_supplement(tid, sys.argv[2])
                    from_supplement = True

                    if lineage is None:
                        lineage = ["Not_in_taxdb"]

                except IndexError:
                    lineage = "Not_in_taxdb"


        #lines.append((name, hit, evalue, lineage))
        if not from_supplement:
            ranks = ncbi.get_rank(lineage)
            ids_at_rank = [[i,r] for [i,r] in ranks.items() if r in rank_of_int]
            for lvl in ids_at_rank:
                lvl[0] = ncbi.get_taxid_translator([lvl[0]])
                lvl[0] = lvl[0][list(lvl[0].keys())[0]]
            ordered = []
            for r in rank_of_int:
                found = False
                for lvl in ids_at_rank:
                    if lvl[1] == r:
                        ordered.append(lvl[0])
                        found = True
                if not found:
                    ordered.append("NA")
            print ('\t'.join([name,hit,evalue,bitscore])+'\t'+'\t'.join(ordered))

        else:
            print ('\t'.join([name,hit,evalue,bitscore])+'\t'+'\t'.join(lineage))

