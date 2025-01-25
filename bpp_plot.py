#!/usr/bin/env python3
"""draw and compare the base pair probabilities of two sequences"""
"""adapted from https://github.com/weiyutang1010/ncrna_design/blob/main/bpp_plot/plot.py"""

import os
import sys
import math
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# from vienna import base_pair_probs
from ensemble import get_bpp_lp

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'


output_dir = "plots/bpps" # default output directory
min_prob = 0.001 # default minimum probability to draw a base pair


def draw_bpps(bpp1, bpp2, pairs_diff, pairs_same, name='beta'):
    seq_len = len(bpp1)
    # Create a new figure for drawing
    figwidth, figheight = max(12, seq_len // 15), 10
    plt.figure(figsize=(figwidth, figheight))

    # Correctly place labels for 5' and 3' ends
    label_dist_x, label_dist_y = seq_len * 0.025, 0
    plt.text(0 - label_dist_x, label_dist_y, r"$5'$", ha='left', va='center', fontsize=12)
    plt.text(seq_len-1 + label_dist_x, label_dist_y, r"$3'$", ha='center', va='center', fontsize=12)

    # Remove axis lines and ticks
    plt.axis('equal')
    plt.axis('off')

    # Place indices on the axis
    # indices_height = seq_len * (-0.010)
    indices_height = -2
    for index in np.linspace(0, seq_len - 1, 10).round().astype(int):
        # plt.text(index, indices_height, str(index + 1), ha='center', va='center', fontsize=13, color='black',)
        plt.text(index, indices_height, fr'${index + 1}$', ha='center', va='center', fontsize=13, color='black',)


    # draw pairs with different probabilities
    def draw_pairs_diff(pairs_diff, bpp, color='blue', scale=1):
        for pair in pairs_diff:
            start_index, end_index, pair_prob = pair[0], pair[1], bpp[pair]
            round_prob = math.ceil(pair_prob * 10) / 10
            # colors = {-1.0: "#FF0000", -0.9: "#FF1A1A", -0.8 : "#FF3333", -0.7: "#FF4D4D", -0.6: "#FF6666", -0.5: "#FF8080", -0.4: "#FF9999", -0.3: "#FFB3B3", -0.2: "#FFCCCC", -0.1: "#FFE6E6", 0.0: "#FFA000", 0.1: "#CCCCFF", 0.2: "#B3B3FF", 0.3: "#9999FF", 0.4: "#8080FF", 0.5: "#6666FF", 0.6: "#4D4DFF", 0.7: "#3333FF", 0.8: "#1A1AFF", 0.9: "#0000FF"}
            # alpha = 1 if not pair_prob else pair_prob
            # equation for ellipses
            # ref: https://en.wikipedia.org/wiki/Ellipse

            center = ((start_index + end_index) / 2, 0)
            width = (end_index - start_index) / 2
            height = (width+3) * 0.55

            x = np.linspace(start_index, end_index, max(400, seq_len * 4))
            y = np.sqrt((1. - ((x - center[0]) ** 2 / (width * width))) * (height * height)) + center[1]

            y *= scale
            if scale == -1:
                y -= 4

            alpha = round_prob
            plt.plot(x, y, alpha=alpha, color=color)

    draw_pairs_diff(pairs_diff, bpp1, color='blue', scale=1)
    draw_pairs_diff(pairs_diff, bpp2, color='blue', scale=-1)
    # draw_pairs_diff(pairs_diff, bpp2, color='blue', scale=-1)

    # draw pairs with similar probabilities
    def draw_pairs_same(pairs_same, bpp, color='green', scale=1):
        for pair in pairs_same:
            start_index, end_index, pair_prob = pair[0], pair[1], bpp[pair]
            global min_prob
            if pair_prob < min_prob:
                continue
            round_prob = math.ceil(pair_prob * 10) / 10
            # colors = {-1.0: "#FF0000", -0.9: "#FF1A1A", -0.8 : "#FF3333", -0.7: "#FF4D4D", -0.6: "#FF6666", -0.5: "#FF8080", -0.4: "#FF9999", -0.3: "#FFB3B3", -0.2: "#FFCCCC", -0.1: "#FFE6E6", 0.0: "#FFA000", 0.1: "#CCCCFF", 0.2: "#B3B3FF", 0.3: "#9999FF", 0.4: "#8080FF", 0.5: "#6666FF", 0.6: "#4D4DFF", 0.7: "#3333FF", 0.8: "#1A1AFF", 0.9: "#0000FF"}
            # color = 'green'
            center = ((start_index + end_index) / 2, 0)
            width = (end_index - start_index) / 2
            height = (width+3) * 0.55

            x = np.linspace(start_index, end_index, max(400, seq_len * 4))
            y = np.sqrt((1. - ((x - center[0]) ** 2 / (width * width))) * (height * height)) + center[1]

            y *= scale

            if scale == -1:
                y -= 4

            alpha = round_prob
            plt.plot(x, y, alpha=alpha, color=color)

    # draw_pairs_same(pairs_same, bpp1, color='orange', scale=1)
    # draw_pairs_same(pairs_same, bpp2, color='orange', scale=-1)
    draw_pairs_same(pairs_same, bpp1, color='blue', scale=1)
    draw_pairs_same(pairs_same, bpp2, color='blue', scale=-1)

    # Draw axis
    x = np.linspace(0, seq_len-1, seq_len)
    y = np.zeros(seq_len)
    plt.plot(x, y, color='black')

    plt.plot(x, y-4, color='black')

    save_path = os.path.join(output_dir, f'{name}.pdf')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    plt.savefig(save_path, bbox_inches='tight')
    print(f"bpps plot saved to {save_path}")

    plt.close()


def compare_bpps(seq1, seq2, name, threshold=0.001):
    # bpp1 = base_pair_probs(seq1)
    # bpp2 = base_pair_probs(seq2)
    print('name:', name)
    print('threshold:', threshold)
    bpp1 = get_bpp_lp(seq1)
    bpp2 = get_bpp_lp(seq2)

    diff = np.abs(bpp1 - bpp2)
    matrix_diff = np.where(diff > threshold)
    matrix_same = np.where(diff <= threshold)
    pairs_diff = [(i, j) for i, j in zip(matrix_diff[0], matrix_diff[1]) if i < j]
    pairs_same = [(i, j) for i, j in zip(matrix_same[0], matrix_same[1]) if i < j]
    assert len(pairs_diff) + len(pairs_same) == len(bpp1) * (len(bpp1) - 1) // 2
    print('diff. pairs size:', len(pairs_diff))
    print('same  pairs size:', len(pairs_same))
    
    draw_bpps(bpp1, bpp2, pairs_diff, pairs_same, name=name)


def main(args):
    df = pd.read_csv(args.input_file)

    # filter rows whose name contains Ours
    df_ours = df[df["name"].str.contains("Ours")] 
    # filter rows whose name contains LinearDesign
    df_ld = df[df["name"].str.contains("LinearDesign")]

    for i in range(len(df_ours)):
        seq1 = df_ld["mRNA"].values[i]
        seq2 = df_ours["mRNA"].values[i]
        name = df_ours["name"].values[i].split("|")[0][1:]
        if name == args.name or args.name == "rna":
            compare_bpps(seq1, seq2, name, args.threshold)
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Draw RNA plots')
    parser.add_argument('--input_file', '-i', default="data/design_results.fasta.csv", type=str, help='Input file')
    parser.add_argument('--name', '-n', default="rna", type=str, help='Name of the sequence to plot')
    parser.add_argument('--out_dir', '-o', default="plots/bpps", type=str, help='Output folder')
    parser.add_argument('--threshold', '-t', type=float, default=0.01, help='Threshold for comparing bpps')
    parser.add_argument('--min_prob', '-m', type=float, default=0.001, help='Minimum probability to draw a base pair')
    parser.add_argument('--online', '-l', action='store_true', help='Use online mode')
    
    args = parser.parse_args()

    output_dir = args.out_dir
    min_prob = args.min_prob

    print("output dir:", output_dir)
    print("min prob:", min_prob)

    if args.online:
        # read input from stdin, two lines per record
        for line1, line2 in zip(sys.stdin, sys.stdin):
            seq1, seq2 = line1.strip(), line2.strip()
            compare_bpps(seq1, seq2, name=args.name, threshold=args.threshold)
    else:
        main(args)