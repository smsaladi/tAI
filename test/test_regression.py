# -*- coding: utf-8 -*-

import os

import pandas as pd
import numpy as np
import Bio.SeqIO

from tAI import tAI

def check_against_ref(fna_file, ref_csv):
    # Read reference features
    calc = pd.read_csv(ref_csv, header=0, names=['old'], squeeze=False)

    my_tai = tAI.from_named_reference('codonR', keep_codonR_err=True)

    with open(fna_file, 'r+') as fh:
        calc['new'] = pd.Series([my_tai.calc(str(rec.seq))
                                    for rec in Bio.SeqIO.parse(fh, 'fasta')])

    assert np.allclose(calc['old'].values, calc['new'].values, atol=0.01)

    return

dir_path = os.path.dirname(os.path.realpath(__file__))

def test_ecolik12():
    basename = os.path.join(dir_path, "data/ecolik12.ffn")
    check_against_ref(basename, basename + ".tAI")
    return


def test_daley_gfp():
    basename = os.path.join(dir_path, "data/Daley_gfp.fna")
    check_against_ref(basename, basename + ".tAI")
    return

test_ecolik12()


"test/data/sacCer3-tRNAs.bed"
