#!/bin/env python

import sys
import os
import regex
from hits.structure.secstruc import read_fasta

_COMP = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "U": "A",
        "N": "N",
        }

def rev_comp(s):
    return "".join([_COMP[c] for c in s.upper()[::-1]])

def fingerprint_search(fingerprint, seq):
    seq = seq.upper()
    _p = regex.compile(r"(%s){s<=1}" % fingerprint.upper())
    ret = []
    for i in range(len(seq)):
        if _p.match(seq[i:]):
            ret.append( (i, _p.match(seq[i:])) )
    return ret

def data_generator(fn):
    with open(fn) as f:
        while True:
            lines = [f.readline() for i in range(4)]
            if not lines[0]:
                break
            yield lines
    return

def processor_builder(sequence_dict, transposon_sequence):
    p1 = regex.compile(r"(%s){s<=1}"%transposon_sequence.upper()[:21])
    p2 = regex.compile(r"(%s){s<=1}"%transposon_sequence.upper()[-21:])

    def processor(lines):
        #read_id = lines[0].strip()
        seq = lines[1].strip()

        # If antisense read
        seq_r = rev_comp(seq)
        if p1.search(seq_r) or p2.search(seq_r):
            seq = seq_r


        # If N-terminus of transposon is in read
        s1 = p1.search(seq)
        if s1:
            seq_v = seq[:s1.span()[0]]
            seq_v_fp = seq_v[-20:] #Use 20 residues to identify ARG and position
            if len(seq_v_fp) == 20:
                for gene in sequence_dict:
                    rseq = sequence_dict[gene]
                    fps1 = fingerprint_search(seq_v_fp, rseq)
                    #If transposon was inserted in parallel
                    if fps1:
                        nt_idx_1st_downstream = fps1[0][0] + len(seq_v_fp) - 5
                        n_rounding = nt_idx_1st_downstream
                        while n_rounding % 3 != 0:
                            n_rounding += 1
                        aa_idx_1st_downstream = (n_rounding / 3) + 1
                        n_extra_nt = n_rounding - nt_idx_1st_downstream
                        return (gene, "P1", aa_idx_1st_downstream, n_extra_nt)
                    #If transposon was inserted in antiparallel
                    else:
                        fps2 = fingerprint_search(rev_comp(seq_v_fp), rseq)
                        if fps2:
                            n_rounding = fps2[0][0]
                            while n_rounding % 3 != 0:
                                n_rounding += 1
                            aa_idx_1st_downstream = (n_rounding / 3) + 1
                            n_extra_nt = n_rounding - fps2[0][0]
                            return (gene, "AP1", aa_idx_1st_downstream, n_extra_nt)
                return "s1"

        # If C-terminus of transposon is in read
        s2 = p2.search(seq)
        if s2:
            seq_v = seq[s2.span()[1]:]
            seq_v_fp = seq_v[:20]
            if len(seq_v_fp) == 20:
                for gene in sequence_dict:
                    rseq = sequence_dict[gene]
                    fps1 = fingerprint_search(seq_v_fp, rseq)
                    #If transposon was inserted in parallel
                    if fps1:
                        n_rounding = fps1[0][0]
                        while n_rounding % 3 != 0:
                            n_rounding += 1
                        aa_idx_1st_downstream = (n_rounding / 3) + 1
                        n_extra_nt = n_rounding - fps1[0][0]
                        return (gene, "P2", aa_idx_1st_downstream, n_extra_nt)
                    #If transposon was inserted in antiparallel
                    else:
                        fps2 = fingerprint_search(rev_comp(seq_v_fp), rseq)
                        if fps2:
                            nt_idx_1st_downstream = fps2[0][0] + len(seq_v_fp) - 5
                            n_rounding = nt_idx_1st_downstream
                            while n_rounding % 3 != 0:
                                n_rounding += 1
                            aa_idx_1st_downstream = (n_rounding / 3) + 1
                            n_extra_nt = n_rounding - nt_idx_1st_downstream
                            return (gene, "AP2", aa_idx_1st_downstream, n_extra_nt)
            return "s2" 
        return "missing"
    return processor


if __name__ == '__main__':
    from hits.sequencing.log_para import parmap, listener_builder
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("fastq", help="Assembled reads fastq")
    parser.add_argument("sequences", help="Functional gene nucleotide sequences")
    parser.add_argument("transposon", help="Inserted design - One sequence!")
    parser.add_argument("outdir", help="Output directory")
    args = parser.parse_args()

    sequence_dict = read_fasta(args.sequences)
    transposon_sequence = list(read_fasta(args.transposon).values())[0]

    listener_process = listener_builder(sequence_dict, args.outdir)

    processor = processor_builder(sequence_dict, transposon_sequence)

    parmap(processor, data_generator(args.fastq), listener_process)

    #Non-parallelized (for debugging)
    '''
    for lines in data_generator(args.fastq):
        data = processor(lines)
        if data:
            print(data)
        else:
            print("None")
        #break
    '''    
