
from Bio import SeqIO
from Bio.pairwise2 import align

ss_dict = {'H': 'H', 'G': 'H', 'I': 'H',
           'B': 'E', 'E': 'E'}

def indexmap(refseq, seq, exact=True):
    if exact:
        alignment = align.globalxs(refseq, seq, -0.5, -0.1)
    else:
        from Bio.SubsMat import MatrixInfo as matlist
        #alignment = align.globalds(refseq, seq, matlist.blosum62, -10, -1)
        alignment = align.localds(refseq, seq, matlist.blosum62, -10, -1)
    #print alignment[0]
    i1, i2 = (0,0)
    out = {}
    for pair in zip(alignment[0][0], alignment[0][1]):
        if pair[0] != '-' and pair[1] != '-':
            out[i2] = i1
            i1+=1
            i2+=1
        elif pair[0] == '-':
            i1+=1
            #print "REFSEQ should not have gap!!"
        elif pair[1] == '-':
            i2+=1
    return out


def _read_dssp(fn):
    pdbseq, ss_elements = [], []
    with open(fn) as f:
        for i, line in enumerate(f):
            if i<25: continue
            aa = line[13]
            if aa == '' or aa == '!': continue
            pdbseq.append(aa)
            if line[16] == ' ': ss = 'C'
            else: ss = line[16]
            ss_elements.append(ss_dict.get(ss,'C'))

    return "".join(pdbseq), ss_elements

def read_fasta(fn):
    data = {}
    with open(fn) as f:
        for record in SeqIO.parse(f, 'fasta'):
            data[record.id] = str(record.seq)
    return data


def get_ss_blocks(dssp_fn, refseq):
    pdbseq, ss_elements = _read_dssp(dssp_fn)
    imap = indexmap(pdbseq, refseq, exact=False)
    
    ss_blocks = []
    prev_ss = None
    idx = 0
    for idx0, this_ss in enumerate(ss_elements):
        if idx0 not in imap: 
            prev_ss = 'D'
            ss_blocks[-1][-1] = idx
            continue
        idx = imap[idx0]+1
        if this_ss != prev_ss:
            ss_blocks.append([this_ss, idx, idx])
            prev_ss = this_ss
        ss_blocks[-1][-1] = idx
    #ss_order = {'H': 1, 'E': 2, 'C': 0}
    #ss_blocks.sort(key=lambda x: _ss_order.get(x[0]))
    return ss_blocks, imap

def read_ddomain(fn, refseq, imap):
    data = []
    index = []
    with open(fn) as f:
        for i, line in enumerate(f):
            ss = line.split()
            if i in imap:
                index.append(imap[i]+1)
                data.append(float(ss[1]))
    return index, data



