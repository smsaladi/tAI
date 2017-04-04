# -*- coding: utf-8 -*-
'''Python adaptation of codonR for tAI calculations

Shyam Saladi (saladi@caltech.edu)
January 2016
'''

# Python 2 and 3 compatibility
from __future__ import print_function

# Only used in main() to read fasta files
# not for tAI calculation
import argparse
import sys

# used in tAI calculation
import re

from pkg_resources import resource_filename

# used in tAI calculation
import pandas as pd
import numpy as np
import scipy.stats.mstats

# Only used in main() to read fasta files
# not for tAI calculation
import Bio.SeqIO

trna_data = None
"""tRNA abundances in references genome. Read upon first invocation.
"""

def trna_table():
    """Read tRNA abundance data from file upon first invocation.
    """
    global trna_data
    if trna_data is None:
        trna_data = pd.read_csv(
            resource_filename(__name__, 'data/codon_abundances.csv'),
                              header=0, index_col=0, comment='#')
    return trna_data

def test_calc_tAI():
    # do tests
    return

def calc_tAI(nucseq, ref_trna='codonR', recalc_weights=False, omit_first=True,
             bacteria=True, optimized_weights=True, keep_codonR_err=True, ):
    """Calculate tAI for the provided nucleotide sequence

    Parameters
    ----------
    nucseq : str

    ref_trna : Optional[str, dict]
        Specifies the reference tRNA abundances used for calculation. `codonR`
        refers to the abundances given with the original package. If specified
        as a `dict`, abundances are taken directly where keys are codons and
        values are counts of cognate tRNAs and values are cached until
        `recalc_weights` = True.

    recalc_weights : Optional[bool]
        Weights for each codon are calculated only once using `calc_weights`
        and then cached/stored in `trna_data`. If `True`, weights will be
        recalculated.

    omit_first : Optional[bool]
        codonR removes the first codon from calculation (whether it's a
        methionine or not). This option is for backwards compatibility.

    bacteria : Optional[bool]
        Passed to `calc_weights`

    optimized_weights : Optional[bool]
        Passed to `calc_weights`

    keep_codonR_err : Optional[bool]
        Passed to `calc_weights`

    Returns
    -------
    float
        tAI value for the nucseq specified with the options given

    Raises
    ------
    None
    """
    # Set up reference tRNA abundance source, i.e. which column
    # or a user supplied dict
    if isinstance(ref_trna, dict):
        lnweights_colname = 'user_ln_weights' + bacteria*'_bact'
        if 'lnweights_colname' not in trna_table().columns or recalc_weights:
            ref_trna_count = pd.Series(ref_trna)
            ref_trna_count.index = \
                ref_trna_count.index.str.upper().replace('U', 'T')
    elif ref_trna in trna_table().columns:
        ref_trna_count = trna_table()[ref_trna]
        lnweights_colname = ref_trna + '_ln_weights' + bacteria*'_bact'
    else:
        raise ValueError('ref_trna not recognized')

    # Now we calculate the relative adaptiveness values for each codon for
    # the reference genome (if not already calculated)
    # Specify a new column to store the ln(weights) given by `lnweights_colname`
    if lnweights_colname not in trna_table().columns or recalc_weights:
        # Only do this once, so store after calculating
        trna_table()[lnweights_colname] = \
            calc_weights(trna_count=ref_trna_count,
                         bacteria=bacteria,
                         optimized_weights=optimized_weights,
                         keep_codonR_err=keep_codonR_err)

    ln_weights = trna_table()[lnweights_colname]
    ln_weights.name = 'ln_weights'

    # We will ignore Methionine codons in our analysis (there is no
    # automatic way to differentiate between 'START' Met-tRNA genes and
    # normal Met-tRNAs in any genome):
    if 'ATG' in ln_weights.index:
        ln_weights.drop('ATG', inplace=True)

    # Check nucseq (just in case)
    nucseq = nucseq.upper().replace('U', 'T')

    if omit_first:
        nucseq = nucseq[3:]

    # Now, count codons:
    # codons = collections.Counter(re.findall('...', nucseq))
    # codon_count = pd.Series(codons, name='codons')
    codon_id, codon_count = np.unique(re.findall('...', nucseq),
                                      return_counts=True)
    codon_count = pd.Series(codon_count, index=codon_id, name='codons')

    # remove non-standard codons from analysis by joining with weights series
    df = pd.concat([ln_weights, codon_count], axis = 1).loc[ln_weights.index]

    # and now we an finally calculate tAI:
    # tai is the weighted average of the weights (weighted by the codon counts)
    codon_prop = df['codons'] / df['codons'].sum()
    return np.exp(np.sum(df['ln_weights'] * codon_prop))

    # This may be faster
    # intermed = np.average(df['ln_weights'], weights=df['codons'], returned=True)
    # return np.exp(intermed[0]/intermed[1])

    # For more details, read the references!
    # [1] dos Reis et al. (2003) Nuc. Acids Res. 31:6976
    # [2] dos Reis et al. (2004) Nuc. Acids Res. 32:5036


def test_calc_weights():
    # do tests
    return

def calc_weights(trna_count, bacteria, optimized_weights, keep_codonR_err):
    """Calculate relative adaptiveness values for each codon

    Parameters
    ----------
    trna_count : pd.Series
        tRNA counts by the codon recognized

    bacteria : Optional[bool]
        If True, consider ATA as coding for isoleucine (with a special
        adaptiveness value).

    optimized_weights : Optional[bool]
        If True, consider optimized (not uniform) weights for each base.

    keep_codonR_err : Optional[bool]
        A transposition in the codonR tRNA reference file resulted in erronous
        weights for TGN codons. If True, this error is kept (for backwards
        compatibility).

    Returns
    -------
    pd.Series
        Relative adaptiveness values for each codon

    Raises
    ------
    None
    """

    if optimized_weights:
        p = {'T': 0.59, 'C': 0.72, 'A': 0.0001, 'G': 0.32}

        # s is what the original script works with
        # p = 1 - s
        # s <- c(0.0, 0.0, 0.0, 0.0, 0.41, 0.28, 0.9999, 0.68)

        isoleucine_p = 1 - 0.89
    else:
        p = {'T': 0.5, 'C': 0.5, 'A': 0.25, 'G': 0.5}
        # s <- c(0, 0, 0, 0, 0.5, 0.5, 0.75, 0.5)
        isoleucine_p = 1 - 0.5

    # obtain absolute adaptiveness values (Ws)

    # add zero's for stop codons since (e.g. 'TGN' may need a 'TGA' abundance)
    trna_count = pd.concat([trna_count,
        pd.Series(0, index=['TGA', 'TAA', 'TAG'])], axis=0)

    # trna file with codonR has a transposition resulting in erronous weights
    # for TGN codons
    if keep_codonR_err:
        trna_count['TGA'] = 1
        trna_count['TGC'] = 1
        trna_count['TGT'] = 0

    # for new weights to be calculated
    weights = pd.Series(0.0, index=trna_count.index)
    # number of tRNAs corresponds to anticodons (not codons)
    for codon in list(weights.index):
        base = codon[:2]
        if codon[2] == 'T':                  # INN -> NNT, NNC, NNA
            weights[codon] = trna_count[codon] + p['T']*trna_count[base+'C']
        elif codon[2] == 'C':                # GNN -> NNT, NNC
            weights[codon] = trna_count[codon] + p['C']*trna_count[base+'T']
        elif codon[2] == 'A':                # TNN -> NNA, NNG
            weights[codon] = trna_count[codon] + p['A']*trna_count[base+'T']
        elif codon[2] == 'G':                # CNN -> NNG
            weights[codon] = trna_count[codon] + p['G']*trna_count[base+'A']
        else:
            raise ValueError('Non-standard codon or notation')

    # Weight calculation in original script
    #   p[1]*tRNA[i] + p[5]*tRNA[i+1],     # INN -> NNT, NNC, NNA
    #   p[2]*tRNA[i+1] + p[6]*tRNA[i],       # GNN -> NNT, NNC
    #   p[3]*tRNA[i+2] + p[7]*tRNA[i],       # TNN -> NNA, NNG
    #   p[4]*tRNA[i+3] + p[8]*tRNA[i+2])     # CNN -> NNG

    # i correspondance based on ordering of input file (i.e. ecolik12.trna)
    # i -> T
    # i+1 -> C
    # i+2 -> A
    # i+3 -> G

    # if bacteria, modify isoleucine ATA codon
    if bacteria:
        weights['ATA'] = isoleucine_p

    # get rid of stop codons and methionine
    weights.drop(['ATG', 'TGA', 'TAA', 'TAG'], inplace=True)

    # get ws
    weights = weights/weights.max()

    # substitute 0-ws by gm
    nonzero_weights = weights[weights != 0]
    geometric_mean = scipy.stats.mstats.gmean(nonzero_weights)
    weights[weights == 0] = geometric_mean

    return np.log(weights)


def main():
    parser = argparse.ArgumentParser(
        description='Calculate tRNA adaptation index')

    parser.add_argument('fna_filename',
        metavar='fna_file',
        type=str,
        default='-',
        help='FASTA-formatted nucleotide file of the Coding Sequences'
             '(stdin by default)')

    args = parser.parse_args()

    if args.fna_filename == '-':
        print("Reading sequences from stdin", file=sys.stderr)
        seq_records = Bio.SeqIO.parse(args.fna_file, "fasta")
    else:
        with open(args.fna_filename, 'r+') as fh:
            seq_records = list(Bio.SeqIO.parse(fh, "fasta"))

    for record in seq_records:
        print(calc_tAI(str(record.seq)))

    return


if __name__ == '__main__':
    main()


## NOT USED FOR tAI CALCULATION
def nc_adj(nc, gc3, a=-6.0, b=34.0, c=1.025):
    """Optimised Nc adjusted function
    a, b, and c have already been optimised
    """
    return a + gc3 + (b/(gc3**2 + (c - gc3)**2)) - nc
