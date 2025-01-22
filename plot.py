#!/usr/bin/env python3
import os
import argparse

import pandas as pd
import matplotlib.pyplot as plt


def plot(file, outfile='entropy_ave'):
    df = pd.read_csv(file)
    # filter rows whose name contains Ours
    df_ours = df[df["name"].str.contains("Ours")] 
    # filter rows whose name contains LinearDesign
    df_ld = df[df["name"].str.contains("LinearDesign")]
    y = df_ours["entropy"]
    x = df_ld["entropy"]

    # plot the data
    _, ax = plt.subplots()
    ax.scatter(x, y, c="blue", alpha=0.5)
    plt.xlabel("LinearDesign")
    plt.ylabel("Ours")

    # make the scale of x and y the same
    ax.set_aspect("equal")

    # draw a diagonal line
    min_val = min(x.min(), y.min())
    max_val = max(x.max(), y.max())
    ax.plot([min_val, max_val], [min_val, max_val], color="red")

    # save the plot to a pdf file
    plt.savefig(f"plots/{outfile}.pdf")
    plt.show()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--file", "-f", type=str, default="data/design_results.fasta.csv", help="path to csv file")
    parser.add_argument("--outfile", "-o", type=str, default="entropy_ave", help="output file name")
    args = parser.parse_args()

    if args.file:
        plot(args.file, args.outfile)
    else:
        print("no file specified")