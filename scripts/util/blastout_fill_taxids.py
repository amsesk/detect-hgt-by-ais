import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-b", "--blastout", required=True, action="store", help="Path to blastout.")
parser.add_argument("-t", "--taxids", required=True, action="store", help="Path to two-column, ;-separated file with accession in column 1 and taxid in column 2")
args = parser.parse_args()

acc_to_tax = {}
with open(args.taxids, 'r') as taxids:
    for line in taxids:
        spl = [x.strip() for x in line.split(";")]
        acc_to_tax[spl[0]] = spl[1]

with open(args.blastout, 'r') as blastout:
    for line in blastout:
        spl = [x.strip() for x in line.split("\t")]

        spl[12] = acc_to_tax[spl[1]]

        print("\t".join(spl))
