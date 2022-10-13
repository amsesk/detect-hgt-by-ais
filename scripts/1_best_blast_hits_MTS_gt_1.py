import sys
from collections import OrderedDict

BITCOL = 11

blastout = open(sys.argv[1], 'r')

current_query = None
current_subject = None

current_query_hits = OrderedDict()

for line in blastout:
    spl = [x.strip() for x in line.split('\t')]
    query = spl[0]
    subject = spl[1]

    # If current_query is None then we are on the first line
    # So set current query equal to the first query
    if current_query is None:
        current_query = query

    # If current_query is equal to query then we're still on the same query
    # So keep adding to the dictionary
    if current_query == query:
        if subject in current_query_hits:
            current_query_hits[subject].append(spl)
        else:
            current_query_hits[subject] = [spl]

    # If current_query is not equal to query, then we've moved onto the next query
    # So pull the best hits to each subject from the dictionary and print them
    else:
        for subject,hit_lines in current_query_hits.items():

            #if len(hit_lines) > 1:
            #    print(len(hit_lines), hit_lines)

    # Modified version of best_blast_hits.py
    # Difference is we're working with a dictionary of subject: [hits to that subject]
    # instead of a file
            best_hit_to_subject = None
            for hit_line in hit_lines:
                if best_hit_to_subject is None:
                    best_hit_to_subject = hit_line

                else:
                    this_bitscore = hit_line[BITCOL]
                    if float(this_bitscore) > float(best_hit_to_subject[BITCOL]):
                        best_hit_to_subject = hit_line

                    else:
                        pass

            print ('\t'.join(best_hit_to_subject))

    # Now that that is done, we can move on to the next query
    # We need to add the current line to the new dictionary here
    # Since we know the dictionary is new, we don't need to test to see
    # if the subject is in the dictionary yet, because it can't be
        current_query_hits = OrderedDict()
        current_query = query

        current_query_hits[subject] = [spl]

