import math
import IsoSpecPy
import numpy as np
from scipy.stats import norm, uniform, gamma
import random
import heapq
import re
from collections import Counter
import numpy.random as rd
from scipy.signal import argrelmax
from .peptides import get_protein_formula
from warnings import warn

from masserstein import Spectrum


class NMRSpectrum(Spectrum):
    def __init__(self, intensities_threshold=0, formula='', threshold=0.001, total_prob=None,
                 charge=1, adduct=None, confs=None, label=None, protons=None, **other):
        if confs is not None:
            confs = np.array(confs)
            if len(confs.shape) != 2 or confs.shape[1] != 2:
                raise ValueError("intensities are not properly shaped")
            confs = confs[confs[:, 1] > intensities_threshold].tolist()
        Spectrum.__init__(self, formula=formula, threshold=threshold, total_prob=total_prob,
                          charge=charge, adduct=adduct, confs=confs, label=label, **other)
        self.protons = protons
