import sys

import mdtraj as mt
from rdkit import Chem


def long2short(code):

    with open("amino-acid.lib") as lines:
        code_mapper = {
            s.split()[2]: s.split()[3] for s in lines if len(s) and s[0] != "#"
        }
    return code_mapper.get(code, "")


def res_seq(pdb):

    p = mt.load_pdb(pdb)
    seq = list(p.topology.residues)
    seq = [str(x)[:3] for x in seq if len(str(x)) >= 3]

    return list(map(long2short, seq))


def lig_smiles(pdb):

    m = Chem.MolFromMol2File(pdb)
    return Chem.MolToSmiles(m)


if __name__ == "__main__":

    p = sys.argv[1]
    o = sys.argv[2]
    s = sys.argv[4]
    lig = sys.argv[3]

    seq = res_seq(p)

    with open(o, "w") as tofile:
        l = f"{p}," + "".join(seq) + "\n"
        tofile.write(l)

    with open(s, "w") as tofile:
        l = f"{lig},{lig_smiles(lig)}"
        tofile.write(l)
