#!/usr/bin/env python3
# call linearpartition to get unpaired probs

# Usage:
# fasta processing: python ensemble.py -f test.fasta -b 0 -d 0
# online mode: echo "GGGCUCGUAGAUCAGCGGUAGAUCGCUUCCUUCGCAAGGAAGAGGCCCUGGGUUCAAAUCCCAGCGAGUCCACCA" | python ensemble.py -b 0 -d 0

import os
import sys
import random
import argparse
import subprocess
from collections import defaultdict


import numpy as np
from scipy.stats import entropy

BASE = 2 # base for entropy calculation


def count_occur(s, char):
    # This function counts the occurrences of 'char' in the string 's'
    count = 0
    for c in s:
        if c == char:
            count += 1
    return count


def call_LP(mRNA, b=0, d=0, keep_bpp=False): # LP (for exact search, use b=0)

    fingerprint = ""
    for nuc in "ACGU":
        count = count_occur(mRNA, nuc)
        fingerprint += f"{nuc}{count}"
    
    hash_value = hash(mRNA)%(10**20)
    rand = random.randint(0, 1000000)
    suffix = f"h{hash_value}_r{rand}_f{fingerprint}_b{b}_d{d}"

    bpp_file = "/tmp/bpp_" + suffix

    # prob = defaultdict(float)
    n = len(mRNA)
    cmd = f"echo {mRNA} | ~/biology/LinearPartition/linearpartition -V -b{b} -d{d} -r {bpp_file}"
    result = subprocess.check_output(cmd, shell=True, text=True) # (((...))) ( -1.20)
    prob = np.identity(n)
    for line in open(bpp_file): # i j prob
        if line.strip() == "":
            break
        i, j, p = line.split()
        i, j, p = int(i)-1, int(j)-1, float(p)
        assert p >= 0
        prob[i, j] = prob[j, i] = p
        prob[i, i] -= p
        prob[j, j] -= p
    prob = np.maximum(prob, 0)
    unp = np.mean(np.diag(prob))
    global BASE
    ent = entropy(prob, base=BASE).mean()
    print("avg unp: %.4f" % unp)
    print("avg ent: %.4f" % ent)
    for i in range(n):
        assert prob[i, i] >= 0, f'''prob[{i}, {i}] = {prob[i, i]}'''
        ent_i = entropy(prob[i])
        # print(f"{i}\t{ent_i:.4f}")
        if ent_i < 0:
            print("negative entropy!")
            print(prob[i])
            print(sum(prob[i]))
    # remove bpp file
    if not keep_bpp and os.path.exists(bpp_file):
        os.remove(bpp_file)
    if keep_bpp:
        print(f"bpp file: {bpp_file}")

    return unp, ent


def batch_call_LP(file, b=0, d=0, keep_bpp=False): # processsing fasta file
    import pandas as pd

    data = []
    for line in open(file):
        if line.startswith(">"): # two lines per record
            name = line.strip()
            continue
        mRNA = line.strip()
        unp, ent = call_LP(mRNA, b, d, keep_bpp)
        data.append((name, mRNA, unp, ent))
    
    df = pd.DataFrame(data, columns=["name", "mRNA", "unpair", "entropy"])
    outfile = f'{file}.csv'
    df.to_csv(outfile, index=False)
    print(f"output to {outfile}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--file", "-f", type=str, default="", help="path to fasta file")
    parser.add_argument("--beam", "-b", type=int, default=0, help="beam size: 0 for exact search")
    parser.add_argument("--dangle", "-d", type=int, default=0, help="dangle mode, 0 or 2")
    parser.add_argument("--bpp", "-p", action="store_true", help="keep bpp files")
    parser.add_argument("--base", "-B", type=float, default=2, help="base for entropy calculation")

    args = parser.parse_args()

    b = args.beam
    d = args.dangle
    print(f"b: {b}, d: {d}")

    BASE = args.base if args.base > 1 else None
    print(f"entropy BASE: {BASE}")

    if args.file:
        batch_call_LP(args.file, b, d, args.bpp)
    else:
        for line in sys.stdin:
            seq = line.strip()
            print(seq)
            call_LP(seq, b, d, args.bpp)