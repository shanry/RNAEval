#!/usr/bin/env python3
import os
import argparse

import pandas as pd
import matplotlib.pyplot as plt


FONT_SIZE_1 = 24
FONT_SIZE_2 = 18


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
    # fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(x, y, c="blue", alpha=0.7)
    plt.xlabel("LinearDesign", fontsize=FONT_SIZE_1)
    plt.ylabel("Ours", fontsize=FONT_SIZE_1)

    # set fontsize of x and y ticks
    plt.xticks(fontsize=FONT_SIZE_2)
    plt.yticks(fontsize=FONT_SIZE_2)

    # make the scale of x and y the same
    ax.set_aspect("equal")

    # tighten layout, remove padding
    plt.tight_layout(pad=0.1)

    # draw a diagonal line
    min_val = min(x.min(), y.min())
    max_val = max(x.max(), y.max())
    ax.plot([min_val, max_val], [min_val, max_val], color="red")

    # Remove excess whitespace with these settings
    # plt.tight_layout(pad=0.0)
    # Adjust the subplot parameters to eliminate margins
    # plt.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1)


    # save the plot to a pdf file
    plt.savefig(f"plots/{outfile}.pdf", bbox_inches='tight', pad_inches=0.05)
    print(f"file saved to plots/{outfile}.pdf")
    # plt.show()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--file", "-f", type=str, default="data/design_results.fasta.csv", help="path to csv file")
    parser.add_argument("--outfile", "-o", type=str, default="entropy_ave", help="output file name")
    args = parser.parse_args()

    if args.file:
        plot(args.file, args.outfile)
    else:
        print("no file specified")