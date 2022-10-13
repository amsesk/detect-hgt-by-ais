import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-b", "--blastout", action = "store", help = "Path to blastout.")
parser.add_argument("--remove", action = "store", help = "Path to newline-separated list of accessions to remove from blastout.")
args = parser.parse_args()

to_remove = []
with open(args.remove) as r:
    for line in r:
        spl = [x.strip() for x in line.split("\t")]
        to_remove.append(spl[0])

with open(args.blastout) as b:
    for line in b:
        spl = [x.strip() for x in line.split("\t")]
        if spl[1] not in to_remove:
            print(line.strip())

