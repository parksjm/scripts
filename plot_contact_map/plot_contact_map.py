# Plot contact maps from a file in CASP RR format
# JMP 4/10/2020
#
# python plot_contact_map.py --con gremlin_50_rr.txt --seq_len 433
# 
# where gremlin_50_rr.txt is something like:
# 
# 147 159 0 8 1.000
# 113 116 0 8 1.000
#  77 161 0 8 1.000
#  28  32 0 8 1.000
# 363 366 0 8 1.000
# 359 400 0 8 1.000
# 181 242 0 8 1.000
# ...

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)


def get_args():
    parser = argparse.ArgumentParser(description='Plot contact maps')
    parser.add_argument("--con", type=str, default="gremlin_50_rr.txt",
                        dest="contactfile", help="contacts in rr format")
    parser.add_argument("--seq_len", type=int, default="433",
                        dest="seq_len", help="number of residues")
    parser.add_argument("--png", type=str, default="test.png",
                        dest="outpng", help="name of output png")

    args = parser.parse_args()

    contactfile = args.contactfile
    seq_len = args.seq_len
    outpng = args.outpng

    return contactfile, seq_len, outpng


def read_couplings(contacts):
    with open(contacts) as sc:
        inum = []
        jnum = []
        prob = []
        for num, line in enumerate(sc, 1):
            # Skip the first line.
            if not num == 1:
                spl = line.split()
                inum.append(int(spl[0]))
                jnum.append(int(spl[1]))
                prob.append(float(spl[4]))

    return inum, jnum, prob 


def plot_contacts3(inum, jnum, probability, seq_len, outpng):
    mtx = np.empty((seq_len, seq_len))
    mtx[:] = np.nan
    for i in range(len(probability)):
        mtx[inum[i], jnum[i]] = probability[i] 
        mtx[jnum[i], inum[i]] = probability[i] 
 
    x, y = np.meshgrid(np.arange(mtx.shape[1]), np.arange(mtx.shape[0]))
    
    fig, ax = plt.subplots(figsize=(3.1, 3.1))

    ax.set_xlim(0,433)
    ax.set_ylim(0,433)

    ax.xaxis.set_major_locator(MultipleLocator(100))
    ax.yaxis.set_major_locator(MultipleLocator(100))

    # Change minor ticks to show every 10. 100/10 = 10
    ax.xaxis.set_minor_locator(AutoMinorLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator(10))

    # Turn grid on for both major and minor ticks and style minor slightly differently
    ax.grid(which='major', color='#CCCCCC', linestyle='--')
    ax.grid(which='minor', color='#CCCCCC', linestyle=':')

    plt.gca().invert_yaxis()
    ax.xaxis.tick_top()

    plt.scatter(x, y, c = mtx, s=2, marker='o', cmap='Blues')
    plt.xlabel("Residue")
    plt.ylabel("Residue")

    plt.savefig(outpng, dpi=300)
    plt.show()
    plt.close()

def main():
    # Get command line arguments
    contactfile, seq_len, outpng = get_args()

    print('contactfile', contactfile)
    print('seq_len', seq_len)
    print('outpng', outpng)

    # Read the coupling scores
    inum, jnum, probability = read_couplings(contactfile)

    # Plot contact map.
    plot_contacts3(inum, jnum, probability, seq_len, outpng)

if __name__ == "__main__":
    main()

