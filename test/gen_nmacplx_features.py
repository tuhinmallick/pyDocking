#!/usr/bin/env python

import numpy as np
import pandas as pd
import mdtraj as mt
import itertools
import sys
from collections import OrderedDict
from mpi4py import MPI
import time


class AtomTypeCounts(object):
    """Featurization of Protein-Ligand Complex based on
    onion-shape distance counts of atom-types.

    Parameters
    ----------
    pdb_fn : str
        The input pdb file name.
    lig_code : str
        The ligand residue name in the input pdb file.

    Attributes
    ----------
    pdb : mdtraj.Trajectory
        The mdtraj.trajectory object containing the pdb.
    receptor_indices : np.ndarray
        The receptor (protein) atom indices in mdtraj.Trajectory
    ligand_indices : np.ndarray
        The ligand (protein) atom indices in mdtraj.Trajectory
    rec_ele : np.ndarray
        The element types of each of the atoms in the receptor
    lig_ele : np.ndarray
        The element types of each of the atoms in the ligand
    lig_code : str
        The ligand residue name in the input pdb file
    pdb_parsed_ : bool
        Whether the pdb file has been parsed.
    distance_computed : bool
        Whether the distances between atoms in receptor and ligand has been computed.
    distance_matrix_ : np.ndarray, shape = [ N1 * N2, ]
        The distances between all atom pairs
        N1 and N2 are the atom numbers in receptor and ligand respectively.
    counts_: np.ndarray, shape = [ N1 * N2, ]
        The contact numbers between all atom pairs
        N1 and N2 are the atom numbers in receptor and ligand respectively.

    """

    def __init__(self, pdb_fn, lig_code):

        self.pdb = mt.load(pdb_fn)

        self.receptor_indices = np.array([])
        self.ligand_indices = np.array([])

        self.rec_ele = np.array([])
        self.lig_ele = np.array([])

        self.lig_code = lig_code

        self.pdb_parsed_ = False
        self.distance_computed_ = False

        self.distance_matrix_ = np.array([])
        self.counts_ = np.array([])

    def parsePDB(self, rec_sele="protein", lig_sele="resname LIG"):

        top = self.pdb.topology

        self.receptor_indices = top.select(rec_sele)
        self.ligand_indices = top.select(lig_sele)

        table, bond = top.to_dataframe()

        self.rec_ele = table['element'][self.receptor_indices]
        self.lig_ele = table['element'][self.ligand_indices]

        self.pdb_parsed_ = True

        return self

    def distance_pairs(self):

        if not self.pdb_parsed_:
            self.parsePDB()

        all_pairs = itertools.product(self.receptor_indices, self.ligand_indices)

        if not self.distance_computed_:
            self.distance_matrix_ = mt.compute_distances(self.pdb, atom_pairs=all_pairs)[0]

        self.distance_computed_ = True

        return self

    def cutoff_count(self, cutoff=0.35):

        self.counts_ = (self.distance_matrix_ <= cutoff) * 1.0

        return self


def generate_features(complex_fn, lig_code, ncutoffs, all_elements):

    # A list of different types of molecules
    #all_elements = ["H", "C", "O", "N", "P", "S", "Br", "Du"]
    #keys = ["_".join(x) for x in list(itertools.product(all_elements, all_elements))]

    cplx = AtomTypeCounts(complex_fn, lig_code)
    cplx.parsePDB(rec_sele="protein", lig_sele="resname %s" % lig_code)

    '''lig = cplx.lig_ele
    rec = cplx.rec_ele

    new_lig, new_rec = [], []
    for e in lig:
        if e not in all_elements:
            new_lig.append("Du")
        else:
            new_lig.append(e)
    for e in rec:
        if e not in all_elements:
            new_rec.append("Du")
        else:
            new_rec.append(e)'''
    # element types of all atoms in the proteins and ligands
    new_lig = [x if x in all_elements else "Du" for x in cplx.lig_ele]
    new_rec = [x if x in all_elements else "Du" for x in cplx.rec_ele]

    rec_lig_element_combines = ["_".join(x) for x in list(itertools.product(new_rec, new_lig))]
    cplx.distance_pairs()

    counts = []
    onion_counts = []

    for i, cutoff in enumerate(ncutoffs):
        cplx.cutoff_count(cutoff)

        if i == 0:
            onion_counts.append(cplx.counts_)
        else:
            onion_counts.append(cplx.counts_ - counts[-1])

        counts.append(cplx.counts_)

    results = []

    for n in range(len(ncutoffs)):
        #count_dict = dict.fromkeys(keys, 0.0)
        d = OrderedDict()
        d = d.fromkeys(keys, 0.0)
        for e_e, c in zip(rec_lig_element_combines, onion_counts[n]):
            d[e_e] += c

        results += d.values()

    return results


if __name__ == "__main__":

    start = time.time()
    print("Start Now ... ")
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # A list of different types of molecules
    all_elements = ["H", "C", "O", "N", "P", "S", "Br", "Du"]
    keys = ["_".join(x) for x in list(itertools.product(all_elements, all_elements))]

    if rank == 0:
        if len(sys.argv) < 3:
            print("Usage: python gen_feature.py input.dat output.csv ")
            sys.exit(0)

        with open(sys.argv[1]) as lines:
            lines = [x for x in lines if ("#" not in x and len(x.split()) >= 1)].copy()
            inputs = [x.split()[0] for x in lines]

        inputs_list = []
        aver_size = int(len(inputs) / size)
        print(size, aver_size)
        for i in range(size-1):
            inputs_list.append(inputs[int(i*aver_size):int((i+1)*aver_size)])
        inputs_list.append(inputs[(size-1)*aver_size:])

        #print(inputs_list)

    else:
        inputs_list = None

    inputs = comm.scatter(inputs_list, root=0)
    #print(rank, inputs)

    out = sys.argv[2]
    n_cutoffs = np.linspace(0.1, 3.1, 60)

    results = []
    ele_pairs =[]
    success = []

    for p in inputs:
        fn = p
        lig_code = "LIG UNK"

        try:
            r = generate_features(fn, lig_code, n_cutoffs, all_elements)
            results.append(r)
            success.append(1.)
            print(rank, fn)

        except:
            #r = results[-1]
            r = list([0., ]*3840) 
            results.append(r)
            #success.append(0.)
            print("Not successful. ", fn)

    df = pd.DataFrame(results)
    try:
        df.index = inputs
    except:
        df.index = np.arange(df.shape[0])

    col_n = []
    #col_n = ['pdbid']
    for i, n in enumerate(keys * len(n_cutoffs)):
        col_n.append(n+"_"+str(i))
#    col_n.append('success')
    df.columns = col_n
#    df['success'] = success
    df.to_csv(str(rank)+"_"+out, sep=",", float_format="%.1f", index=True)

    print(rank, "Complete calculations. ")
    print(time.time() - start)

