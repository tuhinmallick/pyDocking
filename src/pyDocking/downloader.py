# -*- coding: utf-8 -*-
import os
import subprocess as sp
import time
from random import random

import requests
from bs4 import BeautifulSoup
from bs4 import BeautifulSoup as Soup
from pubchempy import get_compounds


class PubChemDownloader(object):
    """Compounds downloader from PubChem.

    Notes
    -----

    Examples
    --------
    >>> # download progesterone and save it into a pdb file
    >>> from pyDocking import downloader
    >>> down = downloader.PubChemDownloader()
    >>> m = down.get_compound("progesterone", type="name")
    >>> smile = down.get_smile(m)
    >>> # now convert SMILE to pdb
    >>> from pyDocking import builder
    >>> b = builder.CompoundBuilder("pdb", "smile")
    >>> b.load_mol(smile)
    >>> b.generate_conformer()
    >>> b.write_mol("progesterone.pdb")

    """

    def __init__(self):
        pass

    def get_compound(self, name, type="name"):
        """Get the compound given its name.

        Parameters
        ----------
        name : str
            The general name of a compound
        type : str, default = name
            The name format.

        Returns
        -------
        results : list
            The list of PubchemPy.compound object
        """

        return get_compounds(name, namespace=type)[0]

    def download_from_list(self, name_list, type="name"):
        """Given a list of compounds, download them

        Parameters
        ----------
        name_list
        type

        Returns
        -------

        """

        results = []
        for n in name_list:
            try:
                r = get_compounds(n, namespace=type)
            except:
                print(f"Ligand {n} not found")
                r = []

            if len(r):
                results.append(r[0])
            else:
                results.append(None)

        return results

    def get_smile(self, compound):
        """Return the smile code of a compound

        Parameters
        ----------
        compound : pcp.Compound object

        Returns
        -------
        smile : str
            The smile code of a molecule

        """

        return compound.isomeric_smiles if compound is not None else ""

    def get_cid(
        self,
        compound,
    ):
        """Return the cid of a compound

        Parameters
        ----------
        compound : pubchempy.compound object

        Returns
        -------
        cid : str
            The cid of a compound.
        """
        return compound.cid if compound is not None else ""


class ZINCDownloader(object):

    def __init__(self):
        pass

    def crawl_smiles(self, url):
        try:
            # url = "http://zinc.docking.org/substance/%s" % zinc_id
            r = requests.get(url)
            soup = Soup(r.text)

            if len(soup.find_all("input")) >= 4:
                return soup.find_all("input")[3]["value"]
            else:
                return ""

        except (ConnectionError, UnicodeDecodeError) as e:
            return ""

    def get_by_id(self, zinc_id):

        if not os.path.exists(zinc_id):

            try:
                cmd = f"wget http://zinc.docking.org/substance/{zinc_id}"
                job = sp.Popen(cmd, shell=True)
                job.communicate()
            except ConnectionError:
                try:
                    url = f"http://zinc.docking.org/substance/{zinc_id}"
                    return self.crawl_smiles(url)

                except (ConnectionError, UnicodeDecodeError) as e:
                    return ""

        if not os.path.exists(zinc_id):
            return ""
        try:
            soup = BeautifulSoup(open(zinc_id))
            if len(soup.find_all("input")) >= 4:
                s = soup.find_all("input")[3]["value"]
            else:
                try:
                    url = f"http://zinc.docking.org/substance/{zinc_id}"
                    s = self.crawl_smiles(url)
                except (ConnectionError, UnicodeDecodeError):
                    s = ""

            return s

        except UnicodeDecodeError:
            try:
                url = f"http://zinc.docking.org/substance/{zinc_id}"

                return self.crawl_smiles(url)
            except (UnicodeDecodeError, ConnectionError) as e:
                return ""

    def get_by_ids(self, zinc_ids, max_sleep=4.0, verbose=True):

        sleep_time = max_sleep * random()
        time.sleep(sleep_time)

        smiles = []
        for zid in zinc_ids:
            try:
                s = self.get_by_id(zid)
                if verbose:
                    print(s, zid)
            except RuntimeError:
                s = ""
            smiles.append(s)

        return smiles
