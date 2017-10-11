# -*- coding: utf-8 -*-
'''Python adaptation of codonR for tAI calculations

Shyam Saladi (saladi@caltech.edu)
January 2016
'''

# Python 2 and 3 compatibility
from __future__ import print_function

# used in tAI calculation
import re

from pkg_resources import resource_filename

# used in tAI calculation
import pandas as pd
import numpy as np
import scipy.stats.mstats

class tAI:
    def __init__(self, trna_count, bacteria=True, optimized_weights=True, keep_codonR_err=False):
        """Initialize a tAI calculator

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
        tAI[object]

        Raises
        ------
        None
        """

        self.weights = self._weights_from_codons(trna_count, bacteria, optimized_weights, keep_codonR_err)
        return


    @staticmethod
    def _weights_from_codons(trna_count, bacteria, optimized_weights, keep_codonR_err):
        """Calculate relative adaptiveness values for each codon.

        See __class__.__init__ for argument decriptions.

        Returns
        -------
        pd.Series
            Relative adaptiveness values for each codon

        Raises
        ------
        ValueError: if a non-ATGC codon is found within the trna_count data
        """

        # Calculate the relative adaptiveness values for each codon for
        # the reference genome (if not already calculated)

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
        weights = pd.Series(0.0, index=trna_count.index, name="ln_weights")
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
        weights.drop(['ATG', 'TGA', 'TAA', 'TAG'], inplace=True, errors='ignore')

        # get ws
        weights = weights/weights.max()

        # substitute 0-ws by geometric mean
        nonzero_weights = weights[~np.isclose(weights, 0)]
        geometric_mean = scipy.stats.mstats.gmean(nonzero_weights)
        weights[weights == 0] = geometric_mean

        return np.log(weights)


    @classmethod
    def from_gtRNAdb(cls, bed_filename, discard_undetermined=True, discard_iMet=True, **kwargs):
        """Initialize a tAI calculator using output files from the GtRNAdb

        See website for more information: http://gtrnadb.ucsc.edu/

        Parameters
        ----------
        bed_filename : str
            tRNA counts by the codon recognized

        discard_undetermined : Optional[bool]
            If True, don't include tRNA's with an undetermined cognate residue
            pairs in the calculation

        discard_iMet : Optional[bool]
            If True, don't include tRNA's that correspond to initiator methionine's
            in the in the calculation

        **kwargs : Optional
            Additional arguments passed to __class__.__init__

        File is expected to be formatted as follows:
            chrI	139151	139254	tRNA-Pro-TGG-1-1	1000	+	139151	139254	0	2	36,36,	0,67,
            chrI	166266	166339	tRNA-Ala-TGC-1-1	1000	+	166266	166339	0	1	73,	0,
            chrIII	149919	149991	tRNA-iMet-CAT-1-1	1000	+	149919	149991	0	1	72,	0,
            chrXII	784353	784453	tRNA-Und-NNN-1-2	1000	+	784353	784453	0	2	36,45,	0,55,

        Returns
        -------
        tAI[object]

        Raises
        ------
        None

        """

        df = pd.read_table(bed_filename, header=None, index_col=False)
        df = df[3].str.split('-', expand=True)[[1, 2]]
        df.columns = ['aa', 'anti']

        if discard_iMet:
            df = df[df['aa'] != 'iMet']

        if discard_undetermined:
            df = df[df['aa'] != 'Und']

        def series_revcomp(ser):
            return (ser.str.replace('A', 't').str.replace('T', 'A').str.replace('t', 'T')
                       .str.replace('G', 'c').str.replace('C', 'G').str.replace('c', 'C')
                       .str[::-1])

        return cls(series_revcomp(df['anti']).value_counts(sort=False), **kwargs)


    @classmethod
    def from_named_reference(cls, name, discard_ATG=True, **kwargs):
        """Initialize a tAI calculator using included abundance data

        Parameters
        ----------
        name : str
            tRNA reference counts for E. coli found within reference file. Either
            `codonR` (the counts provided by the original code) or `lowelab`
            (an updated set of counts from the Lowe Lab's tRNA database, see
            http://gtrnadb.ucsc.edu/genomes/bacteria/Esch_coli_K_12_MG1655/)

        discard_ATG : Optional[bool]
            If True, don't include tRNA's that correspond to methionine's
            in the in the calculation, i.e. as in the original codonR

        **kwargs : Optional
            Additional arguments passed to __class__.__init__

        Returns
        -------
        tAI[object]

        Raises
        ------
        ValueError : `name` not found as a column header in the reference file

        """

        df = pd.read_csv(
            resource_filename(__name__, 'data/codon_abundances.csv'),
            header=0, index_col=0, comment='#')

        if name in df.columns:
            trna_data = df[name]
        else:
            raise ValueError(name, " does not correspond to a reference set provided with this package")

        # Comment from CodonR code:
        # We will ignore Methionine codons in our analysis (there is no
        # automatic way to differentiate between 'START' Met-tRNA genes and
        # normal Met-tRNAs in any genome)
        if discard_ATG and 'ATG' in trna_data.index:
            trna_data.drop('ATG', inplace=True)

        return cls(trna_data, **kwargs)


    def calc(self, seqdata, omit_first=True):
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

        omit_first : Optional[bool]
            codonR removes the first codon from calculation (whether it's a
            methionine or not). This option is for backwards compatibility.

        Returns
        -------
        float
            tAI value for the nucseq specified with the options given

        Raises
        ------
        None
        """
        if isinstance(seqdata, str):
            return self._calc_str(seqdata, omit_first)
        elif isinstance(seqdata, list):
            return self._calc_batch(seqdata, omit_first)
        else:
            raise ValueError(type(seqdata) + " not recognized. Provide `str` or `list(str)`")
        return


    def _calc_str(self, seq, omit_first):
        """Method to calculate tAI for a single sequence
        """
        # Check nucseq (just in case)
        seq = seq.upper().replace('U', 'T')

        if omit_first:
            seq = seq[3:]

        # Now, count codons:
        # codons = collections.Counter(re.findall('...', nucseq))
        # codon_count = pd.Series(codons, name='codons')
        codon_id, codon_count = np.unique(re.findall('...', seq),
                                          return_counts=True)
        codon_count = pd.Series(codon_count, index=codon_id, name='codons')

        # remove non-standard codons from analysis by joining with weights series
        df = pd.concat([self.weights, codon_count], axis=1).loc[self.weights.index]

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


    def _calc_batch(self, seqs, omit_first):
        """Method to calculate tAI scores a list of sequences
        """
        raise NotImplementedError('Batch calculation is not implemented yet, please pass '
                                  'sequences as individual strings')
        return


def main():
    import argparse
    import sys

    import Bio.SeqIO

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

    my_tai = tAI.from_named_reference('codonR')
    print(my_tai.calc([str(r.seq) for r in seq_records]))

    return


if __name__ == '__main__':
    main()
