#!/usr/bin/env python
import os
import sys

import dockml.pdbIO as pdbio


def lig_name_change(lig_in, lig_out, lig_code):

    pio = pdbio.rewritePDB(lig_in)
    with open(lig_out, "w") as tofile:
        with open(lig_in) as lines:
            for s in lines:
                if len(s.split()) and s.split()[0] in ["ATOM", "HETATM"]:
                    nl = pio.resNameChanger(s, lig_code)
                    # n2 = pio.chainIDChanger(nl, "Z")
                    tofile.write(nl)

    return None


def main():

    lig = sys.argv[1]
    out = sys.argv[2]

    with open(lig) as lines:
        lines = [x for x in lines if "LIG" in x]
        if not len(lines):
            lig_name_change(lig, out, "LIG")

        else:
            os.system(f"cp {lig} temp")


main()
