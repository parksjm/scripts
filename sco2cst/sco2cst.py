import sys
import argparse
import math
import numpy as np
import scipy.sparse as sprs
# import pandas as pd
import matplotlib
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import pyplot as plt
from matplotlib import colors as clr
import seaborn as sns; sns.set(color_codes=True)
from pylab import *


def get_args():
    parser = argparse.ArgumentParser(description='Write Rosetta restraints from coevolution scores.')
    parser.add_argument("--msa", type=str, default='PF16998.aln',
                        dest="msa", help="multiple sequence alignment file")
    parser.add_argument("--sco", type=str, default="scores.out",
                        dest="scorefile", help="scores from CCMpred output \
                        after extraction using top_couplings.py")
    parser.add_argument("--sig_cut", type=float, default="0.0",
                        dest="sig_cut", help="probability cutoff for sigmoidal restraints")
    parser.add_argument("--bnd_cut", type=float, default="0.0", dest="bnd_cut",
                        help="probability cutoff for bounded restraints")
    parser.add_argument("--ccm", type=str, default="ccmpred.out", dest="ccm",
                        help="output from CCMpred")

    args = parser.parse_args()

    msa = args.msa
    scorefile = args.scorefile
    sig_cut = args.sig_cut
    bnd_cut = args.bnd_cut
    ccm = args.ccm

    return msa, scorefile, sig_cut, bnd_cut, ccm


def read_alignment(msa):
    try:
        with open(msa) as alignment:
            seq = alignment.readline().rstrip()
    except IOError as err:
        print("I/O error: {0}".format(err))

    return seq

def read_couplings(sco):
    with open(sco) as sc:
        # iaa = []
        inum = []
        jnum = []
        uscore = []
        for num, line in enumerate(sc, 1):
            # Skip the first line.
            if not num == 1:
                spl = line.split()
                inum.append(int(spl[0]))
                jnum.append(int(spl[1]))
                uscore.append(float(spl[2]))

    return inum, jnum, uscore


def scale_scores(x):
    avg_sco = sum(x)/float(len(x))
    min_sco = min(x)
    sscore = []
    for score in x:
        tmpscore = ((score - min_sco)/(avg_sco - min_sco))/2.0 + 0.5
        sscore.append(tmpscore)

    return sscore


def resnum_to_aa(seq, num):
    aa = []
    for resnum in num:
        aa.append(seq[resnum])

    return aa


def get_atomtype(inum, jnum, seq):
    ityp = []
    jtyp = []
    for i, j in zip(inum, jnum):
        if seq[i] == 'G':
            ityp.append('CA')
        else:
            ityp.append('CB')

        if seq[j] == 'G':
            jtyp.append('CA')
        else:
            jtyp.append('CB')

    return ityp, jtyp


def read_cutoffs(iaa, jaa, slope_file):
    slopes = []
    try:
        with open(slope_file) as sc:
            for num, line in enumerate(sc, 1):
                # Skip the first line.
                if not num == 1:
                    splNum = line.split()
                    slopes.append(splNum)
    except IOError as err:
        print("I/O error: {0}".format(err))

    cutoff = []
    sslope = []
    bslope = []

    for i in range(len(iaa)):
        for row in slopes:
            if (row[0] == iaa[i] and row[1] == jaa[i]) \
              or (row[0] == jaa[i] and row[1] == iaa[i]):
                cutoff.append(row[2])
                sslope.append(row[3])
                bslope.append(row[4])

    return cutoff, sslope, bslope


def get_neff(ccm):
    try:
        with open(ccm) as ccmout:
            for line in ccmout:
                if 'Beff' in line:
                    tmp = line.replace('=', ' ').split()
                    beff = tmp[8]
    except IOError as err:
#        print('Failed to open: %s' % err.strerror)
        print("I/O error: {0}".format(err))

    print('neff =', beff)
    return float(beff)


def calc_contact_probs(seq, inum, jnum, nf, sscore):
    N_v_Bx = -0.53
    N_v_By = 5.46
    N_v_Sx = 0.50
    N_v_Sy = 0.58
    N_v_H = 1.0
    N_v_Y = 0.0

    N_v_S = N_v_Sy * nf ** N_v_Sx
    N_v_B = N_v_By * nf ** N_v_Bx

    probs = []
    for i in range(len(inum)):
        p = N_v_H/(1.0 + math.exp(-1.0 * N_v_S * (sscore[i] - N_v_B))) + N_v_Y
        probs.append(p)

    return probs


def plot_contacts(seq, inum, jnum, sscore, probs):
    sns.set_style('white')
    custom_colors = clr.LinearSegmentedColormap.from_list('custom_blue',
                ['#ffffff', '#28e7fd', '#1fb9fc', '#1a98fc', '#1474fb', '#0b28fb'])
    seq_len = int(len(seq))

    mtx = np.zeros((seq_len, seq_len))
    for i in range(len(sscore)):
        mtx[inum[i], jnum[i]] = sscore[i] * probs[i]
        mtx[jnum[i], inum[i]] = sscore[i] * probs[i]

    plt.figure()
    con = sns.heatmap(mtx, cmap=custom_colors, xticklabels=20, yticklabels=20)
    for _, spine in con.spines.items():
        spine.set_visible(True)
    plt.setp(con.xaxis.get_majorticklabels(), rotation=0)
    plt.setp(con.yaxis.get_majorticklabels(), rotation=0)
    con.figure.savefig("contact_map.png", dpi=400)

    print("Wrote contact_map.png")


def plot_contacts2(seq, inum, jnum, sscore, probs):
    seq_len = int(len(seq))

    mtx = np.zeros((seq_len, seq_len))
    # mtx = np.nan((seq_len, seq_len))
    for i in range(len(sscore)):
        mtx[inum[i], jnum[i]] = sscore[i] * probs[i]
        mtx[jnum[i], inum[i]] = sscore[i] * probs[i]
 
    plt.spy(mtx, markersize=1)
    #d=mtx.todense()
    #print(d); sys.exit()
    # https://stackoverflow.com/questions/24013962/how-to-draw-a-matrix-sparsity-pattern-with-color-code-in-python
    # http://scipy-lectures.org/advanced/scipy_sparse/introduction.html
    # d=mtx.todense()
    # plt.imshow(d,interpolation='none', cmap='binary')
    plt.show()


def plot_contacts3(seq, inum, jnum, sscore, probs):
    seq_len = int(len(seq))
    mtx = np.empty((seq_len, seq_len))
    mtx[:] = np.nan
    for i in range(len(sscore)):
        mtx[inum[i], jnum[i]] = sscore[i] * probs[i]
        mtx[jnum[i], inum[i]] = sscore[i] * probs[i]
 
    x, y = np.meshgrid(np.arange(mtx.shape[1]), np.arange(mtx.shape[0]))
    
    plt.axes().set_aspect('equal', 'datalim')
    plt.scatter(x, -y, c = mtx, s=2, marker='s', cmap='copper')
    plt.xlim(0, seq_len)
    plt.ylim(-seq_len, 0)
    plt.show()

def write_restraints(inum, jnum)

    #inum = [x + 1 for x in inum]
    #jnum = [x + 1 for x in jnum]
    #irenum = [x + 1 for x in irenum]
    #jrenum = [x + 1 for x in jrenum]

    try:
        with open('SIG_cst', 'w') as sigOut:
            for i in range(len(irenum)):
                # Scale the score by the probability
                pscore = stuff[4][i] * prob[i]

                if float(stuff[4][i]) >= float(sig_cut):
                    sigOut.write("AtomPair %s %4s %s %4s SCALARWEIGHTEDFUNC %.3f \
                            SUMFUNC 2 SIGMOID %6s %6s CONSTANTFUNC -0.5 #\n"
                            % (stuff[0][i], stuff[1][i], stuff[2][i],
                                stuff[3][i], pscore, stuff[5][i], stuff[6][i]))

def main():
    # Get command line arguments
    msa, scorefile, sig_cut, bnd_cut, ccm = get_args()
    print('msa', msa)
    print('scorefile', scorefile)

    # Read the first sequence in the MSA
    seq = read_alignment(msa)

    # Read the coupling scores
    inum, jnum, uscore = read_couplings(scorefile)

    # Convert the residue numbers to one-letter amino acid symbols
    iaa = resnum_to_aa(seq, inum)
    jaa = resnum_to_aa(seq, jnum)

    # Decide whether to use CB or CA (for Gly only).
    ityp, jtyp = get_atomtype(inum, jnum, seq)

    # Calculate the scaled scores.
    sscore = scale_scores(uscore)

    # Read the cutoff data.
    cutoff, sslope, bslope = read_cutoffs(iaa, jaa, "slopes.txt")

    # Get the effective number of sequences from the ccmpred output.
    neff = get_neff(ccm)
    nf = float(neff)/math.sqrt(float(len(seq)))

    # Calculate the contact probabilities.
    probs = calc_contact_probs(seq, inum, jnum, nf, sscore)

    # Plot contact map.
    plot_contacts3(seq, inum, jnum, sscore, probs)

    # Write out the restraint files
    write_restraints(inum, jnum)

if __name__ == "__main__":
    main()

