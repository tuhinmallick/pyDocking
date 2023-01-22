#!/usr/bin/env python
import os
import sys

from mpi4py import MPI
from pyDocking import builder
from pyDocking import docking
"""
Docking Routine

perform docking with vina with MPI

"""


def do_docking():

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    if rank == 0:
        with open(sys.argv[1]) as lines:
            pdb_codes = [x.strip("\n") for x in lines]
        keep_codes = [
            c
            for c in pdb_codes
            if not os.path.exists(f"{c}/{c}_vinaout.pdbqt")
            and os.path.exists(f"{c}/{c}_protein.pdb.pdbqt")
        ]
        print(keep_codes)
        chunk = int(len(keep_codes) / size)

        input_lists = [
            keep_codes[i * chunk : i * chunk + chunk] for i in range(size - 1)
        ]
        input_lists.append(keep_codes[(size - 1) * chunk:])
    else:
        input_lists = None

    inputs = comm.scatter(input_lists, root=0)

    for c in inputs:

        rec = f"{c}/{c}_protein.pdb.pdbqt"

        # process the ligand
        lig = f"{c}/{c}_ligand.mol2"
        docking.pdb2pdbqt(lig, f"{lig}.pdbqt")
        docking.pdb2pdbqt(lig, f"{lig}.pdb", keep_polarH=False)

        if os.path.exists(rec) and os.path.exists(f"{lig}.pdbqt"):

            try:
                out = f"{c}/{c}_vinaout.pdbqt"
                log = f"{c}/log_vina.log"

                config = f"{c}/vina.config"

                if not os.path.exists(out):
                    rec_prep = docking.ReceptorPrepare(rec)
                    # rec_prep.receptor_addH("H_"+rec)
                    xyz_c = rec_prep.pocket_center(LIG=f"{lig}.pdb")

                    # docking.pdb2pdbqt(rec, "temp.pdbqt")
                    # job = sp.Popen("awk '$1 ~ /ATOM/ {print $0}' temp.pdbqt > %s.pdbqt" % rec, shell=True)
                    # job.communicate()
                    vina = docking.VinaDocking()
                    vina.vina_config(
                        rec,
                        f"{lig}.pdbqt",
                        out,
                        16,
                        8,
                        xyz_c,
                        [20, 20, 20],
                        log,
                        n_modes=20,
                        config=config,
                    )
                    vina.run_docking()

                print("COMPLETE on rank %d: %s" % (rank, c))
            except:
                print("FAIL     on rank %d: %s" % (rank, c))


if __name__ == "__main__":
    do_docking()
