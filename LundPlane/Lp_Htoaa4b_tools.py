import os
import sys
from datetime import datetime
import subprocess
import shlex
import logging
import json
import numpy as np
import math
import awkward as ak
import uproot as uproot
from coffea import hist as coffea_hist
from hist import Hist
from parse import *
import logging
import correctionlib
from coffea import util
from coffea.jetmet_tools import CorrectedJetsFactory, JECStack
from coffea.lookup_tools import extractor
from coffea.nanoevents.methods.nanoaod import GenParticleArray, JetArray
from coffea.nanoevents.methods import candidate, vector
from Lp_Htoaa4b_vars import *

def selectMETFilters(flags_list, era, isMC):
    sFLagDataOrMC = "MC" if isMC else "Data"

    mask_METFilters = np.full(len(flags_list), True, dtype=bool)

    if "goodVertices" in MET_Filters[era][sFLagDataOrMC]:
        mask_METFilters = mask_METFilters & flags_list.goodVertices

    if "globalSuperTightHalo2016Filter" in MET_Filters[era][sFLagDataOrMC]:
        mask_METFilters = mask_METFilters & flags_list.globalSuperTightHalo2016Filter

    if "HBHENoiseFilter" in MET_Filters[era][sFLagDataOrMC]:
        mask_METFilters = mask_METFilters & flags_list.HBHENoiseFilter

    if "HBHENoiseIsoFilter" in MET_Filters[era][sFLagDataOrMC]:
        mask_METFilters = mask_METFilters & flags_list.HBHENoiseIsoFilter

    if "EcalDeadCellTriggerPrimitiveFilter" in MET_Filters[era][sFLagDataOrMC]:
        mask_METFilters = mask_METFilters & flags_list.EcalDeadCellTriggerPrimitiveFilter

    if "BadPFMuonFilter" in MET_Filters[era][sFLagDataOrMC]:
        mask_METFilters = mask_METFilters & flags_list.BadPFMuonFilter

    if "BadPFMuonDzFilter" in MET_Filters[era][sFLagDataOrMC]:
        mask_METFilters = mask_METFilters & flags_list.BadPFMuonDzFilter

    if "hfNoisyHitsFilter" in MET_Filters[era][sFLagDataOrMC]:
        mask_METFilters = mask_METFilters & flags_list.hfNoisyHitsFilter

    if "eeBadScFilter" in MET_Filters[era][sFLagDataOrMC]:
        mask_METFilters = mask_METFilters & flags_list.eeBadScFilter

    if "ecalBadCalibFilter" in MET_Filters[era][sFLagDataOrMC]:
        mask_METFilters = mask_METFilters & flags_list.ecalBadCalibFilter

    return mask_METFilters
