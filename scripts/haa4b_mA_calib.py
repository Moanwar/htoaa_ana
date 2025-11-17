#! /usr/bin/env python
## Calibration study for mass(a) regression in terms of AK8 sub-jets
## Central NanoAOD branches, and custom branches from processor and postprocessor:
## https://cms-nanoaod-integration.web.cern.ch/autoDoc/NanoAODv9/2018UL/doc_SingleMuon_Run2018D-UL2018_MiniAODv2_NanoAODv9-v1.html
## https://gitlab.cern.ch/abrinke1/cmssw/-/blob/PNet_v2_2024_11_22_bkg/PhysicsTools/NanoAOD/python/jets_Hto4b_cff.py
## https://github.com/abrinke1/nanoAOD-tools/blob/nanoPostProc_SS/python/postprocessing/modules/haa4b/objectSelection.py
## https://github.com/abrinke1/nanoAOD-tools/blob/nanoPostProc_SS/python/postprocessing/modules/haa4b/genParticles.py

import os
import sys
import subprocess
import numpy as np
import ROOT as R

R.gROOT.SetBatch(True)  ## Don't display histograms or canvases when drawn

MAX_EVT = -1     ## Maximum number of events to process per MC sample
PRT_EVT = 10000  ## Print every Nth event while processing
DEBUG   = False

## Location of postprocessed input files
IN_DIRS  = '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/YEAR/MC/'
IN_FILES = 'SUSY_GluGluH_01J_HToAATo4B_Pt150_M-mA_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9/r1/PNet_v2_skim.root'


#YEARS  = ['2016APV','2016','2017','2018']
YEARS  = ['2018']
#MASSES = ['15']+[str(15+5*i) for i in range(10)]
#MASSES = ['12']+[str(15+5*i) for i in range(5)]  ## Through 35 GeV
MASSES = ['30']

## Histograms for mass(a) distributions
hst = {}
## Thresholds for Lo/Med/Hi
thr = {}
## Count events passing different selection cuts
count = {}
## Color codes for histograms
color = {}
color['Med_Med'] = R.kBlack
color['Hi_Hi']   = R.kBlue
color['Lo_Lo']   = R.kRed
color['Hi_Med']  = R.kCyan-3
color['Med_Hi']  = R.kViolet
color['Lo_Med']  = R.kMagenta
color['Med_Lo']  = R.kGreen+1
color['Hi_Lo']   = R.kOrange+8
color['Lo_Hi']   = R.kTeal+7

if not os.path.exists('plots/png/haa4b_mA_calib'):
    os.makedirs('plots/png/haa4b_mA_calib')

for year in YEARS:
    print('\n\n*** Beginning to look at year %s ***\n' % year)
    hst[year] = {}
    thr[year] = {}
    count[year] = {}
    for mA in MASSES:
        for var in ['jms','jes','jrs','sjms','sjes','sjrs']:
            ## 2D jet/sub-jet mass and energy and mass/pT scale histograms
            hst[year]['mA_%s_%s2vs1' % (mA, var)] = R.TH2D('h_2D_mA_%s_%s_%s' % (mA, var, year),
                                                           'm(a) = %s, %s (%s)' % (mA, var, year),
                                                           60, 0, 3.0, 60, 0, 3.0)
            for idx in ['1','2']:
                ## 1D jet/sub-jet mass and energy and mass/pT scale histograms
                hst[year]['mA_%s_%s%s' % (mA, var, idx)] = R.TH1D('h_1D_mA_%s_%s%s_%s' % (mA, var, idx, year),
                                                                  'm(a) = %s, %s%s (%s)' % (mA, var, idx, year),
                                                                  60, 0, 3.0)
                thr[year]['mA_%s_%s%s' % (mA, var, idx)] = [-1,-1]
            for js1 in ['Hi','Med','Lo']:
                for js2 in ['Hi','Med','Lo']:
                    ## Categorize by overlapping AK4 jet ("j") or AK8 soft-drop sub-jet ("sj")
                    ## mass scale ("ms") or energy scale ("es") relative to GEN mass(a)
                    cat = '%s_%s1_%s_%s2_%s' % (mA, var, js1, var, js2)
                    hst[year]['mA_'+cat] = R.TH1D('h_1D_mA_'+cat+'_'+year,
                                                  'm(a) = %s, %s1 %s %s2 %s (%s)' % (mA, var, js1, var, js2, year),
                                                  min(int(mA)*5, 70), 5, min(float(mA)+max(5, int(float(mA)/3)), 75))
                    hst[year]['mA_'+cat].SetLineColor(color['%s_%s' % (js1, js2)])
    ## End loop: for js1 in ['Hi','Med','Lo']

    ## Combine input files into a single "chain"
    in_file_names = [IN_DIRS.replace('YEAR',year)+IN_FILES.replace('mA',mA) for mA in MASSES]
    chains = {}
    chains['Events'] = 0
    for i in range(len(in_file_names)):
        print('Adding file %s' % in_file_names[i])
        for key in chains.keys():
            if i == 0: chains[key] = R.TChain(key)
            chains[key].Add( in_file_names[i] )

    ## Loop through events, select, and count
    ch = chains['Events']  ## Shortcut expression
    nEntries = ch.GetEntries()
    print('\nEntering loop over %d events\n' % (nEntries))

    for key in ['None','iH','X4b','4B','iA','2j','2sj']:
        count[year][key] = 0
    for iEvt in range(nEntries):

        if iEvt >= MAX_EVT and MAX_EVT > 0: break
        if (iEvt % PRT_EVT) == 0: print('Looking at event #%d / %d' % (iEvt, nEntries))

        ch.GetEntry(iEvt)
        count[year]['None'] += 1

        ## Good Haa4b candidate AK8 must exist
        iH = ch.Haa4b_iFatH
        if iH < 0: continue
        if ch.FatJet_pt[iH] < 250.0: continue
        count[year]['iH'] += 1
        ## Higgs AK8 must pass X4b WP60
        if ch.FatJet_PNet_X4b_v2a_Haa4b_score[iH] + \
           ch.FatJet_PNet_X4b_v2b_Haa4b_score[iH] < 2.0*0.93: continue  ## Higgs signal must pass WP60
        count[year]['X4b'] += 1
        ## Higgs AK8 must contain 4 GEN-level b-hadrons and quarks
        if ch.FatJet_nBHadrons[iH] < 4: continue
        #if ch.Haa4b_FatH_nBQuarks  < 4: continue
        count[year]['4B'] += 1
        ## Event must have good GEN "a" bosons
        iA1 = ch.GEN_a1_idx
        iA2 = ch.GEN_a2_idx
        if iA1 < 0 and iA2 < 0: continue
        count[year]['iA'] += 1
        mAs  = [int(ch.GenPart_mass[iA1]), int(ch.GenPart_mass[iA2])]
        ptAs = [ch.GenPart_pt[iA1], ch.GenPart_pt[iA2]]
        assert(mAs[0] == mAs[1])
        assert(str(mAs[0]) in MASSES)
        mA = mAs[0]

        ## Look for two highest-pT AK4 jets overlapping Higgs AK8
        jmass = []
        jpt   = []
        for ij in range(ch.nJet):
            if ch.Jet_Haa4b_ovlp_iFat[ij] != iH: continue
            if ch.Jet_Haa4b_presel[ij]    != 1:  continue
            jmass.append(ch.Jet_mass[ij])
            jpt  .append(ch.Jet_pt[ij])
            if len(jmass) == 2: break

        ## Look for Higgs AK8 sub-jets
        sjmass = []
        sjpt   = []
        if ch.FatJet_subJetIdx1[iH] >= 0:
            sj1 = ch.FatJet_subJetIdx1[iH]
            sjmass.append(ch.SubJet_mass[sj1])
            sjpt  .append(ch.SubJet_pt[sj1])
        if ch.FatJet_subJetIdx2[iH] >= 0:
            sj2 = ch.FatJet_subJetIdx2[iH]
            sjmass.append(ch.SubJet_mass[sj2])
            sjpt  .append(ch.SubJet_pt[sj2])

        #########################################################################
        ## Fill histograms based on relation between AK4 jet RECO and GEN mass/pT
        if len(jmass) == 2:
            count[year]['2j'] += 1
            
            ## Fill base jms/jes/jrs histograms
            jsr = {'jms':[], 'jes':[], 'jrs':[]}
            for xj in [0,1]:
                jsr['jms'].append(jmass[xj] / mAs[xj])
                jsr['jes'].append( jpt[xj] / ptAs[xj])
                jsr['jrs'].append(jsr['jms'][xj] / jsr['jes'][xj])
                for var in ['jms','jes','jrs']:
                    key = 'mA_%s_%s%d' % (mA, var, xj+1)
                    hst[year][key].Fill(min(jsr[var][xj], 2.999))
                    if xj == 1:
                        key2 = 'mA_%s_%s2vs1' % (mA, var)
                        hst[year][key2].Fill(min(jsr[var][0], 2.999), min(jsr[var][1], 2.999))
                    ## Update threshold every 1000 events
                    if (int(hst[year][key].Integral()) % 1000) == 0:
                        p = np.array([0.25,0.75])
                        q = np.array([0.,0.])
                        hst[year][key].GetQuantiles(2, q, p)
                        thr[year][key] = [q[0], q[1]]
                        if DEBUG: print('Updated %s thresholds to [%.4f, %.4f]' % (key, q[0], q[1]))

            ## Fill Hi/Lo/Med plots
            for var in ['jms','jes','jrs']:
                ## Fill plots, but only if thresholds are set (after first 1000 events)
                thrL = [thr[year]['mA_%s_%s%d' % (mA, var, xj)][0] for xj in [1,2]]
                thrH = [thr[year]['mA_%s_%s%d' % (mA, var, xj)][1] for xj in [1,2]]
                if thrL[0] < 0 or thrL[1] < 0 or thrH[0] < 0 or thrH[1] < 0: continue
                js1 = 'Hi' if jsr[var][0] > thrH[0] else 'Lo' if jsr[var][0] < thrL[0] else 'Med'
                js2 = 'Hi' if jsr[var][1] > thrH[1] else 'Lo' if jsr[var][1] < thrL[1] else 'Med'
                cat = '%d_%s1_%s_%s2_%s' % (mA, var, js1, var, js2)
                hst[year]['mA_'+cat].Fill(ch.FatJet_PNet_34massAa[iH])
        ## End conditional: if len(jmass) == 2


        ########################################################################
        ## Fill histograms based on relation between subjet RECO and GEN mass/pT
        if len(sjmass) == 2:
            count[year]['2sj'] += 1

            ## Fill base jms/jes/jrs histograms
            sjsr = {'sjms':[], 'sjes':[], 'sjrs':[]}
            for xj in [0,1]:
                sjsr['sjms'].append(sjmass[xj] / mAs[xj])
                sjsr['sjes'].append( sjpt[xj] / ptAs[xj])
                sjsr['sjrs'].append(sjsr['sjms'][xj] / sjsr['sjes'][xj])
                for var in ['sjms','sjes','sjrs']:
                    key = 'mA_%s_%s%d' % (mA, var, xj+1)
                    hst[year][key].Fill(min(sjsr[var][xj], 2.999))
                    if xj == 1:
                        key2 = 'mA_%s_%s2vs1' % (mA, var)
                        hst[year][key2].Fill(min(sjsr[var][0], 2.999), min(sjsr[var][1], 2.999))
                    ## Update threshold every 1000 events
                    if (int(hst[year][key].Integral()) % 1000) == 0:
                        p = np.array([0.25,0.75])
                        q = np.array([0.,0.])
                        hst[year][key].GetQuantiles(2, q, p)
                        thr[year][key] = [q[0], q[1]]
                        if DEBUG: print('Updated %s thresholds to [%.4f, %.4f]' % (key, q[0], q[1]))

            ## Fill Hi/Lo/Med plots
            for var in ['sjms','sjes','sjrs']:
                ## Fill plots, but only if thresholds are set (after first 1000 events)
                thrL = [thr[year]['mA_%s_%s%d' % (mA, var, xj)][0] for xj in [1,2]]
                thrH = [thr[year]['mA_%s_%s%d' % (mA, var, xj)][1] for xj in [1,2]]
                if thrL[0] < 0 or thrL[1] < 0 or thrH[0] < 0 or thrH[1] < 0: continue
                sjs1 = 'Hi' if sjsr[var][0] > thrH[0] else 'Lo' if sjsr[var][0] < thrL[0] else 'Med'
                sjs2 = 'Hi' if sjsr[var][1] > thrH[1] else 'Lo' if sjsr[var][1] < thrL[1] else 'Med'
                cat = '%d_%s1_%s_%s2_%s' % (mA, var, sjs1, var, sjs2)
                hst[year]['mA_'+cat].Fill(ch.FatJet_PNet_34massAa[iH])

        ## End conditional: if len(sjmass) == 2


    ## End loop: for iEvt in range(nEntries)
## End loop: for year in YEARS

print('\n*** Finished looking at all years! ***\n\n')

for year in YEARS:
    nTot = count[year]['None']
    nPrev = nTot
    print('\n%s had %d total events' % (year, nTot))
    for key in ['iH','X4b','4B','iA','2j','2sj']:
        nX = count[year][key]
        print('  * %s cut: %d (%.2f%% lost, %.2f%% total)' % (key, nX, 100*(nPrev-nX)/nPrev, 100*nX/nTot))
        nPrev = nX if key != '2j' else nPrev
## End loop: for year in YEARS

print('\nWriting out to plots/haa4b_mA_calib.root')
outf = R.TFile('plots/haa4b_mA_calib.root', 'recreate')
outf.cd()
cans = {}
maxs = {}
for year in YEARS:
    cans[year] = {}
    maxs[year] = {}
    for key in hst[year].keys():
        hist = hst[year][key]
        keyA = key.replace('Hi','All').replace('Med','All').replace('Lo','All')
        hist.SetLineWidth(2)
        hist.Write()
        can = R.TCanvas(hist.GetName())
        can.cd()
        if key.endswith('2vs1'):
            R.gStyle.SetOptStat(0) ## Don't display stat boxes
            hist.Draw('colz')
        else:
            R.gStyle.SetOptStat(1110) ## Display rms, mean, and number of entries
            hist.Draw('histe')
        can.SaveAs('plots/png/haa4b_mA_calib/'+hist.GetName()+'.png')
        del can

        R.gStyle.SetOptStat(0) ## Don't display stat boxes
        if keyA.count('All') == 2 and not keyA in cans[year].keys():
            maxs[year][keyA] = hist.GetMaximum()
            cans[year][keyA] = R.TCanvas('%s_%s' % (keyA, year))
            cans[year][keyA].cd()
            hist.Draw('histe')
        elif keyA.count('All') == 2:
            maxs[year][keyA] = max(maxs[year][keyA], hist.GetMaximum())
            cans[year][keyA].cd()
            hist.Draw('histesame')
    ## End loop: for key in hst[year].keys()

    R.gStyle.SetOptStat(0) ## Don't display stat boxes
    for keyA in cans[year].keys():
        cans[year][keyA].cd()
        leg = R.TLegend(0.12, 0.42, 0.48, 0.88)
        for js1 in ['Hi','Med','Lo']:
            for js2 in ['Hi','Med','Lo']:
                key = keyA.replace('1_All','1_'+js1).replace('2_All','2_'+js2)
                hist = hst[year][key]
                hist.SetTitle(hist.GetTitle().replace('1 '+js1,'1 All').replace('2 '+js2,'2 All'))
                hist.GetYaxis().SetRangeUser(0, 1.2*maxs[year][keyA])
                leg.AddEntry(hist, '%.2f #pm %.2f (%3s-%3s)' % (hist.GetMean(), hist.GetStdDev(), js1, js2), 'l')
        leg.Draw()
        cans[year][keyA].SaveAs('plots/png/haa4b_mA_calib/%s_%s.png' % (keyA, year))
        del leg

## End loop: for year in YEARS
        
outf.Write()
outf.Close()

print('\n*** All done!!! ***\n\n')
