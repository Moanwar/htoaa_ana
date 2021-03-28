#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 28 14:56:36 2021

@author: si_sutantawibul1
"""
import numpy as np
from EfficiencyInfoB import EfficiencyInfoB


class EfficiencyInfoC(EfficiencyInfoB):
    def __init__(self, demDf, name, var, nbins, histrange):
        super().__init__(demDf, name, var, nbins, histrange)
        self.plotdir = 'JetHTTrigEff/plots/C/'

    def getNumDf(self, df):
        trig = np.logical_and(np.logical_or(df.L1_DoubleJet112er2p3_dEta_Max1p6,
                                            df.L1_DoubleJet150er2p5),
                              df.HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71)
        ret = df[trig]
        return ret