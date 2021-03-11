#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 20:01:45 2021

@author: si_sutantawibul1
"""
import numpy as np
import pandas as pd
import sys
import pickle
from scipy.special import ndtri
from scipy.stats import norm
import matplotlib.pyplot as plt
import dataVsMC_DataManager as DM

class EfficiencyInfo(object) :
    def __init__(self, demDf, name):
        if not (isinstance(demDf, pd.DataFrame)):
            print('init error: input needs to be dataframes')
            sys.exit()

        self.hasWeights = ('final_weights' in demDf)
        self.nbins = 27
        self.range = (150, 1500)
        self.name = name
        self.demDf = demDf
        self.numDf = self.getNumDf(demDf)
        self.dem, self.demEdge = self.hist(demDf)
        self.num, self.numEdge = self.hist(self.numDf)
        self.quot = np.divide(self.num, self.dem, where=self.dem!=0)
        self.upErr, self.lowErr = self.computeError()



        #self.plot()


    ## makes the pt into histograms
    def hist(self, df):
        if self.hasWeights:
            return np.histogram(df['FatJet_pt'], weights=df['final_weights'],
                                bins = self.nbins, range=self.range)
        else:
            return np.histogram(df['FatJet_pt'], bins=self.nbins, range=self.range)

    def getNumDf(self, df):
        return df.loc[(df.L1_SingleJet180==True) & (df.HLT_AK8PFJet500==True)]


    def getBinCenter(self, arr):
        arrCen = list()
        for i in range(len(arr)-1):
            arrCen.append((arr[i+1]+arr[i])/2)
        return arrCen

    ## compute normal error for 1 bin
    ## returns np array of error. Code found at:
    ## https://root.cern.ch/doc/master/TEfficiency_8cxx_source.html#l02744
    def normalError(self, total, passed, weights):
        if 0 == total:
            return 0,0
        if self.hasWeights:
            totalErr = np.sqrt((weights*weights).sum())
        #else:
        #    totalErr =  np.sqrt(len(total))

        level = 0.68
        alpha = (1-level)/2
        avgWgt = np.power(totalErr,2)/total
        jitter = avgWgt/total
        average = passed/total
        sigma = np.sqrt((avgWgt*(average+jitter)*(1+jitter-average))/total)
        delta = norm.ppf(1-alpha,0,sigma)
        #delta = -sigma*ndtri(1-alpha)

        upper = min(delta,1-average)#(average + delta) if ((average + delta) < 1) else 1
        lower = min(delta,average)#(average - delta) if ((average - delta) > 0) else 0

        return upper, lower

    ## compute error for whole histogram
    def computeError(self, ):
        edges = self.demEdge
        demdf = self.demDf
        upperError = list()
        lowerError = list()

        for i in range(len(edges)-1):
            if self.hasWeights:
                wgts = demdf.final_weights.loc[(demdf.FatJet_pt >= edges[i])
                                            & (demdf.FatJet_pt < edges[i+1])]
            else:
                wgts = 0

            tmpUp, tmpLow = self.normalError(self.dem[i], self.num[i], wgts)
            upperError.append(tmpUp)
            lowerError.append(tmpLow)


        return np.array(upperError), np.array(lowerError)

    def plotpt(self, ):
        edge = self.getBinCenter(self.demEdge)

        fig, ax = plt.subplots(figsize=(10,6))
        ax.grid()
        ax.set_ylim([-0.05,1.05])
        ax.set_title(self.name)

        ax.errorbar(edge, self.quot, yerr=(self.lowErr, self.upErr),
                    linestyle='None',fmt='ok', capsize=5)
        xerr = np.ones((2, len(self.quot)))*26
        ax.errorbar(edge, self.quot, xerr=xerr, linestyle='None', fmt='k')

        plotdir = 'JetHTTrigEff/plots/'
        plt.savefig(f'{plotdir}{self.name}.png')