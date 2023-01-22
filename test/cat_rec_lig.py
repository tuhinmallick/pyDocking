#!/usr/bin/env python
import os
import subprocess as sp
import sys

import dockml.pdbIO as pdbio
from pyDocking import builder
from rdkit import Chem


def convert_lig(lig_in, lig_out):

    try:
        mol = Chem.MolFromMol2File(lig_in, removeHs=False)
        Chem.MolToPDBFile(mol, lig_out)
    except:
        print("Switch to open babel. ")
        builder.babel_converter(lig_in, lig_out)

    return None


def cat_rec_lig(rec, lig, out):
    job = sp.Popen("cat %s | awk '$1 ~ /ATOM/ {print $0}' > temp_rec" % rec,
                   shell=True)
    job.communicate()
    job = sp.Popen(
        "cat temp_rec %s | awk '$1 ~ /ATOM/ || $1 ~ /HETATM/ {print $0}' | awk '$4 != /HOH/ {print $0}' > %s"
        % (lig, out),
        shell=True,
    )
    job.communicate()


def lig_name_change(lig_in, lig_out, lig_code):

    pio = pdbio.rewritePDB(lig_in)
    with open(lig_out, "w") as tofile:
        with open(lig_in) as lines:
            for s in lines:
                if len(s.split()) and s.split()[0] in ["ATOM", "HETATM"]:
                    nl = pio.resNameChanger(s, lig_code)
                    n2 = pio.chainIDChanger(nl, "Z")
                    tofile.write(n2)

    return None


def main():

    inputs = [
        x.split()[0] for x in open(sys.argv[1]).readlines() if "#" not in x
    ]
    print(inputs)
    for p in inputs:

        rec = os.path.join(p, f"{p}_protein.pdb")
        lig = os.path.join(p, f"{p}_ligand.mol2")
        try:
            convert_lig(lig, f"t1_{p}.pdb")
            lig_name_change(f"t1_{p}.pdb", f"t2_{p}.pdb", "LIG")
            cat_rec_lig(rec, f"t2_{p}.pdb", os.path.join(p, f"{p}_cplx.pdb"))
        except:
            print("Not successful : ", p)


#        print(p)

main()
