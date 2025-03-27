#!/usr/bin/env python3

"""
This script is used to plot the energy and probabilities of suboptimal structures
"""

import os

import sys

import subprocess

import numpy as np
from matplotlib import pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

FONT_SIZE_1 = 27
FONT_SIZE_2 = 18

FONT_SIZE_LEGEND = 29


def get_free_energy(sequence):
    # Construct the command
    cmd = f"echo {sequence} | ~/biology/LinearPartition/linearpartition -V -d 0 -b 0"
    
    # Run the command and capture the output
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    # Check if the command was successful
    if result.returncode != 0:
        raise RuntimeError(f"Command failed with error: {result.stderr}")
    
    # Combine stdout and stderr (in case the output is in stderr)
    full_output = result.stdout + result.stderr
    
    # Split the output into lines
    output_lines = full_output.strip().split('\n')
    # print(output_lines)
    
    # Extract the free energy value
    for line in output_lines:
        # print(line)
        if "Free Energy of Ensemble:" in line:
            free_energy = float(line.split(":")[1].strip().split()[0])
            return free_energy
    
    # If no free energy line is found
    raise ValueError("Free Energy of Ensemble not found in the output")


def partition(energy_list, cutoff=10):
    # sum of (exp(-1/0.6163207755 * (energy )))
    partition = 0
    for energy in energy_list[:cutoff]:
        partition += np.exp(-1/0.6163207755 * energy)
    return partition


def get_probs(energy_list, cutoff=10):
    partition_value = partition(energy_list, cutoff=cutoff)
    probs = []
    for energy in energy_list:
        prob = np.exp(-1/0.6163207755 * energy) / partition_value
        probs.append(prob)
    return probs


def get_full_probs(seq, energy_list):
    energy_ensemble = get_free_energy(seq)
    energy_list = [x - energy_ensemble for x in energy_list]
    partition_value = np.exp(-1/0.6163207755 * 0)
    probs = []
    for energy in energy_list:
        prob = np.exp(-1/0.6163207755 * energy) / partition_value
        probs.append(prob)
    return probs


def plot_energy(data1, data2, name = 'energy'):
    # Create a figure and axis
    fig, ax = plt.subplots()

    # Plot the data
    ax.plot(data1, label='data1', color='blue')
    ax.plot(data2, label='data2', color='red')

    # set y limit, upper bound
    plt.ylim(min(data1 + data2), -357.1)


    # set title
    # plt.title(f'Energy landscape of suboptimal structures')

    # Add a legend
    ax.legend(['EnsembleDesign', 'LinearDesign'], fontsize=FONT_SIZE_LEGEND)

    # label y unit as kcal/mol
    # plt.ylabel('kcal/mol')
    plt.ylabel(r'$\Delta G^\circ ~\text{(kcal/mol)}$', fontsize=FONT_SIZE_1)
    plt.xlabel('rank of structures in the ensemble', fontsize=FONT_SIZE_1)

    # set fontsize of x and y ticks
    plt.xticks(fontsize=FONT_SIZE_2)
    plt.yticks(fontsize=FONT_SIZE_2)
    
    # tighten layout, remove padding
    plt.tight_layout(pad=0.1)


    # Show the plot
    # plt.show()

    # Save the plot
    fig.savefig(f'plots/energy/{name}.pdf')
    # plt.close()
    print(f'file saved to plots/energy/{name}.pdf')


def plot_probs(data1, data2, cutoff=50, name = 'prob'):

    # minus min. energy to avoid underflow
    data1 = [x - min(data1) for x in data1]
    data2 = [x - min(data2) for x in data2]

    x_values_1 = range(1, len(data1) + 1)
    x_values_2 = range(1, len(data2) + 1)

    # Calculate the probabilities
    probs1 = get_probs(data1, cutoff=cutoff)
    probs2 = get_probs(data2, cutoff=cutoff)

    # Create a figure and axis
    fig, ax = plt.subplots()

    # Plot the data
    ax.plot(x_values_1, probs1, label='probs1', color='blue')
    ax.plot(x_values_2, probs2, label='probs2', color='red')

    # Add a legend
    ax.legend(['EnsembleDesign', 'LinearDesign'])

    # add a title
    plt.title(f'Boltzmann probabilities of {cutoff} suboptimal structures')
    
    # label y unit as p(y|x)
    plt.ylabel( 'normalized ' + r'$p(\boldsymbol{y} \mid \boldsymbol{x})$ ' + f'in the partial ensemble')

    # plt.xlabel('rank of suboptimal structures (' + r'$\boldsymbol{y}$' + ') in ensemble ' + r'$\mathcal{Y}(\boldsymbol{x})$')
    plt.xlabel('rank of structures in the ensemble')

    # Show the plot
    # plt.show()

    # Save the plot
    fig.savefig(f'plots/prob/{name}-{cutoff}.pdf')
    print(f'file saved to plots/prob/{name}-{cutoff}.pdf')
    plt.close()


def plot_full_probs(data1, data2, seq1, seq2, name = 'prob'):
    # Calculate the probabilities
    probs1 = get_full_probs(seq1, data2)
    probs2 = get_full_probs(seq2, data2)

    # return

    # Create a figure and axis
    fig, ax = plt.subplots()

    x_values = range(1, len(data1) + 1)

    # Plot the data
    ax.plot(x_values, probs1, label='probs1', color='blue')
    ax.plot(x_values, probs2, label='probs2', color='red')

    # # tick x axis from 1
    # plt.xticks(np.arange(1, len(probs1)+1, 1.0))

    # set title
    plt.title(f'')

    # Add a legend
    ax.legend(['EnsembleDesign', 'LinearDesign'], fontsize=FONT_SIZE_LEGEND)

    # add a title
    # plt.title(f'Boltzmann probabilities of suboptimal structures') 

    # label y unit as p(y|x)
     # label y unit as p(y|x)
    plt.ylabel(r'$p(\boldsymbol{y} \mid \boldsymbol{x})$ ' + 'in the ensemble', fontsize=FONT_SIZE_1)

    # plt.xlabel('rank of suboptimal structures (' + r'$\boldsymbol{y}$' + ') in ensemble ' + r'$\mathcal{Y}(\boldsymbol{x})$')
    plt.xlabel('rank of structures in the ensemble', fontsize=FONT_SIZE_1)

    # set fontsize of x and y ticks
    plt.xticks(fontsize=FONT_SIZE_2)
    plt.yticks(fontsize=FONT_SIZE_2)

    # tighten layout, remove padding
    plt.tight_layout(pad=0.1)

    # Show the plot
    # plt.show()

    # Save the plot
    fig.savefig(f'plots/exact_prob/{name}.pdf')
    print(f'file saved to plots/exact_prob/{name}.pdf')
    plt.close()


def test(file1, file2):

    # Read in the file
    with open(sys.argv[1], 'r') as f:
        data1 = f.readlines()

    with open(sys.argv[2], 'r') as f:
        data2 = f.readlines()

    # Convert the data to floats
    data1 = [float(x.strip()) for x in data1]
    data2 = [float(x.strip()) for x in data2]

    # Create a figure and axis
    fig, ax = plt.subplots()

    # Plot the data
    ax.plot(data1, label='data1', color='blue')
    ax.plot(data2, label='data2', color='red')

    # Add a legend
    ax.legend(['Ours', 'LinearDesign'])

    # label y unit as kcal/mol
    plt.ylabel(r'$\Delta G ~\text{(kcal/mol)}$')

    # Show the plot
    # plt.show()

    # Save the plot
    fig.savefig('energy.pdf')

    cut = 10

    # Calculate the probabilities
    probs1 = get_probs(data1, cutoff=cut)
    probs2 = get_probs(data2, cutoff=cut)

    # Create a figure and axis
    fig, ax = plt.subplots()

    # Plot the data
    ax.plot(probs1, label='probs1', color='blue')
    ax.plot(probs2, label='probs2', color='red')

    # Add a legend
    ax.legend(['Ours', 'LinearDesign'])

    # add a title
    plt.title(f'cutoff at {cut}')
    
    # label y unit as p(y|x)
    plt.ylabel('p(y|x)')

    # Show the plot
    # plt.show()

    # Save the plot
    fig.savefig('probs.pdf')


def get_ids(dir='data/subopt'):
    ids = set()
    for file in os.listdir(dir):
        ids.add(file.split('.')[0])
    return ids


def main():
    ids = get_ids()
    for id in ids: # ['Q9P2M1']
        if id == 'SPIKE':
            continue
        if id != 'P0DMU9': # P0DMU9 Q13794
            continue
        try:
            print(f'Processing {id}')
            file1 = f'data/subopt/{id}.ours.dG.d0.subopt.3'
            file2 = f'data/subopt/{id}.LD.dG.d0.subopt.3'
            file1_seq = f'data/subopt/{id}.ours'
            file2_seq = f'data/subopt/{id}.LD'

            max_cut = 200

            # Read in the file
            with open(file1, 'r') as f:
                data1 = f.readlines()[:max_cut]

            with open(file2, 'r') as f:
                data2 = f.readlines()[:max_cut]

            with open(file1_seq, 'r') as f:
                seq1 = f.read().strip()

            with open(file2_seq, 'r') as f:
                seq2 = f.read().strip()

            # Convert the data to floats
            data1 = [float(x.strip()) for x in data1[:max_cut]]
            data2 = [float(x.strip()) for x in data2[:max_cut]]

            # show size of data
            print(f'len(data1) = {len(data1)}')
            print(f'len(data2) = {len(data2)}')

            # max_cut = min(len(data1), len(data2), max_cut)

            plot_energy(data1, data2, name=id)
            if max_cut >= max(len(data1), len(data2)):
                plot_probs(data1, data2, cutoff=max_cut, name=id)
            plot_full_probs(data1, data2, seq1=seq1, seq2=seq2, name=id)
        except:
            print(f'Error processing {id}')
            continue

        # break
        
    return


if __name__ == '__main__':

    # file1 = sys.argv[1] # each line is an energy value
    # file2 = sys.argv[2]

    # test(file1, file2)

    main()
