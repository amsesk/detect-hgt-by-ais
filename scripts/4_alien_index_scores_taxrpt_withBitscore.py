import sys
import numpy as np

blastout = sys.argv[1]

alien_index_cutoff = 20
alien_id = "Eukaryota"
domestic_id = "Bacteria"

query_hits = {}
with open(blastout) as b:
    for line in b:
        spl = [x.strip() for x in line.split('\t')]
        if spl[0] in query_hits:
            query_hits[spl[0]].append(spl)
        else:
            query_hits[spl[0]] = [spl]

for k,v in query_hits.items():
    alien = [x for x in v if x[4] == alien_id]
    domestic = [x for x in v if x[4] == domestic_id]
    alien_e = [float(x[2]) for x in alien]
    domestic_e = [float(x[2]) for x in domestic]

    best_aliens = []
    for i,ae in enumerate(alien_e):
        if ae == min(alien_e):
            best_aliens.append(i)

    if len(alien_e) == 0:
        alien_best = 1
    else:
        alien_best = min(alien_e)

    if len(domestic_e) == 0:
        domestic_best = 1
    else:
        domestic_best = min(domestic_e)

    alien_index = np.log((domestic_best)+1e-200) - np.log((alien_best)+1e-200)
    if alien_index >= alien_index_cutoff:
        if len(best_aliens) > 1:
            print (f"[WARNING] More than one best alien hit for: {k}")
        print('\t'.join([k, str(v[best_aliens[0]][1]), str(alien_best), str(domestic_best), str(alien_index)]))
    '''
    else:
        if domestic_best is None:
            print(k, alien_best, domestic_best, np.nan)

        else:
            pass
            #print(k, alien_best, domestic_best, np.nan)
    '''
