import numpy as np
from collections import defaultdict

def read_data(bdir, gene):
    data = defaultdict(dict)
    for condition in ['T', 'A', 'R']:
        tmp_p = []
        tmp_ap = []
        for offset in [0,1,2]:
            x, y_p = np.loadtxt(
                    bdir+"/%s_counts/%s_P_%d.count"%(condition, gene, offset),
                    dtype=int, unpack=True)
            tmp_p.append(y_p[0:-1])
            _, y_ap = np.loadtxt(
                    bdir+"/%s_counts/%s_AP_%d.count"%(condition, gene, offset),
                    dtype=int, unpack=True)
            tmp_ap.append(y_ap[0:-1])
        data[condition]['Parallel'] = np.array(tmp_p)
        data[condition]['Antiparallel'] = np.array(tmp_ap)
    return data
