import os
import sys

import numpy as np
import pandas as pd
from pyDocking import builder
from rdkit import Chem
from rdkit.Chem.Fingerprints import FingerprintMols


def get_top_features(fn, converter):

    m = converter(fn)

    try:
        s = FingerprintMols.FingerprintMol(m).ToBitString()
        f = [float(x) for x in list(s)]
    except:
        f = np.zeros(2048)

    return f


if __name__ == "__main__":

    Mol = builder.Molecule()

    inp = sys.argv[1]
    out = sys.argv[2]

    if not os.path.exists(inp):
        print(
            "python gen_top_feat.py input_pdb_code.dat output_top_features.csv"
        )
        sys.exit(0)

    with open(inp) as lines:
        fn_list = [x.split()[0] for x in lines]

    features = []
    for fn in fn_list:
        lig = f"{fn}/{fn}_ligand.mol2"

        f = get_top_features(lig, Mol.converter_["mol2"])
        features.append(f)

    df = pd.DataFrame(features)
    df.columns = [f"FP{str(x)}" for x in range(df.shape[1])]
    df.to_csv(out, header=True, index=False, sep=",", float_format="%.1f")
