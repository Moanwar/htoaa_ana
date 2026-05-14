#! /usr/bin/env python
## Create ROOT TTrees with VBF events passing X4b selections
## Rewritten in the Coffea framework (coffea >= 2023.x / NanoEvents style)

import os
import sys
import json
import argparse
import urllib.request
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import awkward as ak
import numpy as np
import uproot

from coffea import processor
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
from coffea.lumi_tools import LumiMask

try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False

# ──────────────────────────────────────────────────────────────────────────────
# Constants
# ──────────────────────────────────────────────────────────────────────────────

DEFAULT_YEAR = "2018"

JSON_FILES = {
    "2016preVFP":  "list_samples_Pre2016.json",
    "2016postVFP": "list_samples_Post2016.json",
    "2017":        "list_samples_2017.json",
    "2018":        "list_samples_2018.json",
}

GOLDEN_JSON_URLS = {
    "2016preVFP":  "https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt",
    "2016postVFP": "https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt",
    "2017":        "https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt",
    "2018":        "https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt",
}

VBF_CATEGORIES  = ["VBFjjLoPtLo", "VBFjjLoPtHi", "VBFjjHiPtLo", "VBFjjHiPtHi"]
X4B_SELECTIONS  = ["X4bSR", "X4bSB", "X4bSB_only"]

# Trigger lists per year / primary dataset
TRIGGERS = {
    "2018": {
        "JetHT": [
            "AK8PFJet500", "PFHT1050", "PFJet500",
            "AK8PFHT800_TrimMass50", "AK8PFJet400_TrimMass30",
        ],
        "BTagCSV": [
            "AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4",
            "DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71",
            "QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2",
            "QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1",
            "PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5",
        ],
    },
    "2017": {
        "JetHT": [
            "PFJet500",
            "PFHT380_SixPFJet32_DoublePFBTagCSV_2p2",
            "PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2",
            "PFHT430_SixPFJet40_PFBTagCSV_1p5",
            "PFHT1050",
            "AK8PFHT750_TrimMass50",
            "AK8PFHT800_TrimMass50",
            "AK8PFJet500",
            "AK8PFJet360_TrimMass30",
            "AK8PFJet380_TrimMass30",
            "AK8PFJet400_TrimMass30",
        ],
        "BTagCSV": [
            "DoublePFJets100MaxDeta1p6_DoubleCaloBTagCSV_p33",
            "PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0",
        ],
    },
    "2016": {

        "JetHT": [
            "PFJet450", "DiCentralPFJet430",
            "PFHT650_WideJetMJJ900DEtaJJ1p5", "PFHT750_4JetPt50",
            "PFHT800", "PFHT900", "AK8PFJet360_TrimMass30",
            "AK8PFJet450", "AK8PFHT650_TrimR0p1PT0p03Mass50",
            "AK8PFHT700_TrimR0p1PT0p03Mass50",
            "AK8DiPFJet250_200_TrimMass30_BTagCSV_p20",
            "AK8DiPFJet280_200_TrimMass30_BTagCSV_p20",
            "PFHT400_SixJet30_DoubleBTagCSV_p056",
            "PFHT450_SixJet40_BTagCSV_p056",
            "AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20",
        ],
        "BTagCSV": [
            "DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6",
            "DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160",
            "DoubleJetsC112_DoubleBTagCSV_p014_DoublePFJetsC112MaxDeta1p6",
            "DoubleJetsC112_DoubleBTagCSV_p026_DoublePFJetsC172",
            "DoubleJet90_Double30_TripleBTagCSV_p08",
            "QuadJet45_TripleBTagCSV_p087",
            "QuadPFJet_BTagCSV_p016_VBF_Mqq460",
            "QuadPFJet_BTagCSV_p016_VBF_Mqq500",
            "QuadPFJet_BTagCSV_p016_p11_VBF_Mqq200",
            "QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240",
        ],
    },
}

# ──────────────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────────────

def extract_year_from_dataset(name):
    if "2016" in name:
        return "2016preVFP" if ("HIPM" in name or "Pre" in name or "pre" in name) else "2016postVFP"
    if "2017" in name:
        return "2017"
    if "2018" in name:
        return "2018"
    print(f"Warning: cannot determine year from '{name}'. Using {DEFAULT_YEAR}")
    return DEFAULT_YEAR


def download_golden_json(era):
    url = GOLDEN_JSON_URLS[era]
    print(f"Downloading golden JSON: {url}")
    try:
        with urllib.request.urlopen(url) as r:
            return r.read().decode()
    except Exception as exc:
        print(f"Error downloading golden JSON for {era}: {exc}")
        sys.exit(1)


def get_file_list(dataset_name, year):
    json_file = JSON_FILES.get(year)
    if not json_file or not os.path.exists(json_file):
        print(f"Error: JSON file {json_file} not found")
        return []
    with open(json_file) as f:
        filesets = json.load(f)
    if dataset_name not in filesets:
        print(f"Error: '{dataset_name}' not in {json_file}")
        print(f"Available: {list(filesets.keys())}")
        return []
    return filesets[dataset_name]


def extract_primary_dataset(paths):
    for p in paths:
        for name in ("JetHT", "BTagCSV", "SingleMuon", "SingleElectron", "MET"):
            if name in p:
                return name
    return None


def _trig_key(year):
    """Return the TRIGGERS sub-dict key for a given era string."""
    for k in ("2018", "2017", "2016"):
        if k in year:
            return k
    return None


# ──────────────────────────────────────────────────────────────────────────────
# Coffea processor
# ──────────────────────────────────────────────────────────────────────────────

class VBFx4bProcessor(processor.ProcessorABC):
    """
    Coffea processor that replicates the original event selection and writes
    per-category TTree output via uproot.

    Parameters
    ----------
    year            : era string, e.g. "2018" / "2016preVFP"
    primary_dataset : "JetHT" | "BTagCSV" | None
    lumi_mask       : coffea LumiMask object (or None to skip)
    output_file     : path to the output ROOT file that will be written
    """

    def __init__(self, year, primary_dataset, lumi_mask, output_file):
        self.year            = year
        self.primary_dataset = primary_dataset
        self.lumi_mask       = lumi_mask
        self.output_file     = output_file
        self._trig_key       = _trig_key(year)

        # bTag WP per year
        # DeepJet medium WP: 2018=0.2783, 2017=0.3040, 2016preVFP=0.2598, 2016postVFP=0.2489
        self.btagWPM = {"2018": 0.2783, "2017": 0.3040}.get(
            self._trig_key, 0.2598 if "preVFP" in year else 0.2489
        )

        # Accumulators – simple dict of lists that we'll flush to ROOT
        self._out = {
            f"{c}_{s}": {k: [] for k in (
                "run", "luminosityBlock", "event",
                "fatjet_pt", "fatjet_eta", "fatjet_phi", "fatjet_mass",
                "vbfjet1_pt", "vbfjet1_eta", "vbfjet1_phi", "vbfjet1_mass",
                "vbfjet2_pt", "vbfjet2_eta", "vbfjet2_phi", "vbfjet2_mass",
                "vbfjj_deta", "vbfjj_mass", "x4b_score",
                "pnet_massH_v2b", "pnet_34massAa",
                "trigger_fired", "trigger_priority",
            )}
            for c in VBF_CATEGORIES for s in X4B_SELECTIONS
        }

        self.cnt = {}  # filled dynamically by _cf() in process()

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _safe_field(collection, field, default_val=-999.0):
        """Return collection field, or a constant-filled array if absent."""
        if hasattr(collection, field):
            return getattr(collection, field)
        return ak.full_like(collection.pt, default_val)

    @staticmethod
    def _pt_field(collection):
        return collection.pt_nom

    @staticmethod
    def _mass_field(collection):
        return collection.mass_nom

    @staticmethod
    def _msoftdrop_field(fj):
        return fj.msoftdrop_nom

    @staticmethod
    def _dijet_mass(j1_pt, j1_eta, j1_phi, j1_m, j2_pt, j2_eta, j2_phi, j2_m):
        """Proper invariant mass of two jets via 4-vector sum."""
        px1 = j1_pt * np.cos(j1_phi)
        py1 = j1_pt * np.sin(j1_phi)
        pz1 = j1_pt * np.sinh(j1_eta)
        E1  = np.sqrt(px1**2 + py1**2 + pz1**2 + j1_m**2)

        px2 = j2_pt * np.cos(j2_phi)
        py2 = j2_pt * np.sin(j2_phi)
        pz2 = j2_pt * np.sinh(j2_eta)
        E2  = np.sqrt(px2**2 + py2**2 + pz2**2 + j2_m**2)

        m2 = (E1+E2)**2 - (px1+px2)**2 - (py1+py2)**2 - (pz1+pz2)**2
        return np.sqrt(ak.where(m2 > 0, m2, 0.0))

    @staticmethod
    def _delta_r(eta1, phi1, eta2, phi2):
        deta = eta1 - eta2
        dphi = (phi1 - phi2 + np.pi) % (2 * np.pi) - np.pi
        return np.sqrt(deta**2 + dphi**2)

    # ------------------------------------------------------------------
    # Trigger helpers (return per-event boolean arrays)
    # ------------------------------------------------------------------

    def _hlt_branches(self, events):
        """
        Return the set of HLT branch names actually present in this chunk.
        coffea NanoEvents does NOT support hasattr() on HLT fields — we must
        inspect the awkward array fields directly.
        """
        hlt = events.HLT
        # ak.fields() returns the list of fields for a record array
        return set(ak.fields(hlt))

    def _fired_any(self, events, names, available=None):
        """OR over a list of HLT branch names; missing branches are skipped."""
        if available is None:
            available = self._hlt_branches(events)
        result = ak.zeros_like(events.run, dtype=bool)
        hlt = events.HLT
        matched = [n for n in names if n in available]
        for name in matched:
            result = result | getattr(hlt, name)
        return result, matched

    def _trigger_mask(self, events):
        """
        Returns (pass_mask, priority_code_array).
        priority_code: 1 = JetHT, 2 = BTagCSV, 0 = unknown/2018-single-dataset
        """
        tkey = self._trig_key
        pd   = self.primary_dataset

        if tkey is None:
            print("  WARNING: unknown trig_key, no triggers will fire!")
            return ak.zeros_like(events.run, dtype=bool), ak.zeros_like(events.run, dtype=np.int32)

        available = self._hlt_branches(events)

        # First chunk: report which of our trigger names are found in the file
        if not hasattr(self, "_trig_branches_reported"):
            self._trig_branches_reported = True
            wanted_jht  = set(TRIGGERS[tkey]["JetHT"])
            wanted_btag = set(TRIGGERS[tkey]["BTagCSV"])
            found_jht   = wanted_jht  & available
            found_btag  = wanted_btag & available
            miss_jht    = wanted_jht  - available
            miss_btag   = wanted_btag - available
            print(f"  [trig branches] JetHT  found={sorted(found_jht)}  missing={sorted(miss_jht)}")
            print(f"  [trig branches] BTagCSV found={sorted(found_btag)}  missing={sorted(miss_btag)}")
            if not found_jht and not found_btag:
                # Print a sample of what IS in HLT to help diagnose naming
                sample = sorted(available)[:30]
                print(f"  [trig branches] NONE matched! Sample of HLT fields in file: {sample}")

        fired_jht,  matched_jht  = self._fired_any(events, TRIGGERS[tkey]["JetHT"],  available)
        fired_btag, matched_btag = self._fired_any(events, TRIGGERS[tkey]["BTagCSV"], available)

        if "2018" in self.year:
            mask = fired_jht | fired_btag
            prio = ak.where(fired_jht, np.int32(1), np.int32(2))
            prio = ak.where(mask, prio, np.int32(0))
            return mask, prio

        # 2016 / 2017 – priority-based deduplication
        if pd == "JetHT":
            mask = fired_jht #| fired_btag
            prio = ak.where(fired_jht, np.int32(1), np.int32(2))
            prio = ak.where(mask, prio, np.int32(0))
        elif pd == "BTagCSV":
            mask = fired_btag & (~fired_jht)
            prio = ak.where(mask, np.int32(2), np.int32(0))
        #else:
        #    mask = fired_jht | fired_btag
        #    prio = ak.zeros_like(events.run, dtype=np.int32)

        return mask, prio

    # ------------------------------------------------------------------
    # Selection helpers
    # ------------------------------------------------------------------

    def _noise_flags(self, events):
        """
        MET filters following:
        https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2

        2018 / 2017 Data+MC:
            goodVertices, globalSuperTightHalo2016Filter, HBHENoiseFilter,
            HBHENoiseIsoFilter, EcalDeadCellTriggerPrimitiveFilter,
            BadPFMuonFilter, BadPFMuonDzFilter, eeBadScFilter,
            ecalBadCalibFilter
            (hfNoisyHitsFilter intentionally excluded)

        2016 (preVFP + postVFP) Data+MC:
            same as above EXCEPT no ecalBadCalibFilter,
            PLUS hfNoisyHitsFilter
        """
        F    = events.Flag
        year = self.year

        # Common filters for all eras
        base = (
            F.goodVertices
            & F.globalSuperTightHalo2016Filter
            & F.HBHENoiseFilter
            & F.HBHENoiseIsoFilter
            & F.EcalDeadCellTriggerPrimitiveFilter
            & F.BadPFMuonFilter
            & F.BadPFMuonDzFilter
            & F.eeBadScFilter
        )

        if "2016" in year:
            # 2016: add hfNoisyHitsFilter, no ecalBadCalibFilter
            if hasattr(F, "hfNoisyHitsFilter"):
                base = base & F.hfNoisyHitsFilter
        else:
            # 2017 / 2018: add ecalBadCalibFilter, no hfNoisyHitsFilter
            base = base & F.ecalBadCalibFilter

        return base
    def _muon_veto(self, events):
        mu = events.Muon
        #pt_thr = 26.0 + (3.0 if "2017" in self.year else 0.0)
        pt_thr = 26.0 
        good = (
            (mu.pt > pt_thr)
            & (abs(mu.eta) < 2.4)
            & (mu.mediumPromptId | ((mu.pt > 200) & (mu.highPtId > 0)))  # ← Change mediumId to mediumPromptId
            & (abs(mu.dxy) < 0.02)
            & (abs(mu.dz)  < 0.10)
            & (mu.miniPFRelIso_all < 0.10)
        )
        return ~ak.any(good, axis=1)


    def _electron_veto(self, events):
        ele = events.Electron
        #pt_thr = 30.0 + (5.0 if "2016" not in self.year else 0.0)
        pt_thr = 35.0 
        in_crack = (abs(ele.eta) > 1.44) & (abs(ele.eta) < 1.57)
        
        # Explicitly convert to bool and handle missing values
        wp90_ok = (ele.mvaFall17V2Iso_WP90 >= 1) if hasattr(ele, 'mvaFall17V2Iso_WP90') else False
        wp80_ok = (ele.mvaFall17V2Iso_WP80 >= 1) if hasattr(ele, 'mvaFall17V2Iso_WP80') else False
        heep_ok = (ele.cutBased_HEEP >= 1) if hasattr(ele, 'cutBased_HEEP') else False
        
        good = (
            (ele.pt > pt_thr)
            & (abs(ele.eta) < 2.5) & ((abs(ele.eta) < 1.44) | (abs(ele.eta) > 1.57))
            #& (~in_crack)
            & wp90_ok
            & (wp80_ok | ((ele.pt > 35) & heep_ok))
            & (abs(ele.dxy) < 0.02)
            & (abs(ele.dz)  < 0.10)
        )
        return ~ak.any(good, axis=1)

    def _select_fatjet(self, events):
        """
        Mirrors the reference fat-jet loop exactly:

          for iFat:
              if pt <= 250 / |eta| >= 2.4 / jetId != 6 / msoftdrop <= 20: skip
              check top tagger              <- on base-passing jets (no pt>300 cut)
              check WZ tagger               <- on base-passing jets
              if Xbb <= 0.75: skip          <- only for Higgs candidate
              compute X4b score

        Returns four progressive per-event masks + best-candidate kinematics:
          pass_candH  — base (pt>250, |eta|<2.4, jetId==6, msoftdrop>20) + Xbb>0.75
          pass_leadPt — same as pass_candH (pt>250 already inside base, kept for
                        cutflow labelling only — will equal pass_candH)
          pass_notop  — pass_candH AND no OTHER base-passing fat jet passes top tagger
          pass_nowz   — pass_notop AND no OTHER base-passing fat jet passes WZ tagger
        """
        fj   = events.FatJet
        year = self.year

        fj_pt   = self._pt_field(fj)
        fj_mass = self._mass_field(fj)
        fj_sd   = self._msoftdrop_field(fj)

        # Base quality — INCLUDES pt > 250, matching reference loop
        base = (
            (fj_pt > 250)
            & (abs(fj.eta) < 2.4)
            & (fj.jetId == 6)
            & (fj_sd > 20)
        )


        # Higgs candidate: base + Xbb > 0.75
        xbb_ok    = fj.particleNetMD_XbbvsQCD > 0.75
        cand_mask = base & xbb_ok

        # X4b score: average of Haa4b v2a+v2b (matching reference code)
        if hasattr(fj, "PNet_X4b_v2a_Haa4b_score") and hasattr(fj, "PNet_X4b_v2b_Haa4b_score"):
            x4b_scores = 0.5 * (fj.PNet_X4b_v2a_Haa4b_score + fj.PNet_X4b_v2b_Haa4b_score)
        elif hasattr(fj, "X4b_score"):
            x4b_scores = fj.X4b_score
        else:
            x4b_scores = ak.full_like(fj_pt, -999.0)

        score_for_max = ak.where(cand_mask, x4b_scores, -999.0)

        # Best candidate per event (highest X4b score among cand_mask jets)
        has_any_cand  = ak.any(cand_mask, axis=1)
        n_fj          = ak.num(fj_pt)
        safe_score    = ak.where(n_fj > 0, score_for_max, ak.Array([[0.0]] * len(events)))
        best_idx_flat = ak.argmax(safe_score, axis=1)

        def _gather(arr):
            return arr[ak.local_index(arr) == best_idx_flat[:, None]]

        best_x4b   = ak.flatten(_gather(score_for_max), axis=None)
        best_pt    = ak.flatten(_gather(fj_pt),          axis=None)
        best_eta   = ak.flatten(_gather(fj.eta),         axis=None)
        best_phi   = ak.flatten(_gather(fj.phi),         axis=None)
        best_mass  = ak.flatten(_gather(fj_mass),        axis=None)

        best_massH_v2b = ak.flatten(_gather(fj.PNet_massH_v2b), axis=None)                          if hasattr(fj, "PNet_massH_v2b") else ak.full_like(best_pt, -999.0)
        best_34massAa  = ak.flatten(_gather(fj.PNet_34massAa),  axis=None)                          if hasattr(fj, "PNet_34massAa")  else ak.full_like(best_pt, -999.0)

        # pass_candH: base (incl. pt>250) + Xbb + valid X4b score
        pass_candH  = has_any_cand & (best_x4b > -999.0)

        # pass_leadPt: same as pass_candH since pt>250 is already in base.
        # Kept as a separate mask so the cutflow can still print "leadingFatJetPt".
        pass_leadPt = pass_candH   # pt>250 already enforced in base

        # dR from best Higgs candidate (needed for both top and WZ veto pools)
        dR_from_best = self._delta_r(fj.eta, fj.phi,
                                     best_eta[:, None], best_phi[:, None])

        # nonHto4bFatJet for TOP veto:
        #   selFatJets with dR > 0.05, pt > 300, 50 < msoftdrop < 220
        #   (mirrors: nonHto4bFatTJet = selFatJets[dR>0.05 & pt>300 & 50<sd<220])
        is_top_pool = (
            base                        # pt>250, |eta|<2.4, jetId==6, sd>20
            & (dR_from_best > 0.05)
            & (fj_pt > 300)
            & (fj_sd > 50)
            & (fj_sd < 220)
        )

        # Top veto: leading TvsQCD score in top pool > threshold
        top_thresh = {"2018": 0.970, "2017": 0.970}.get(
            self._trig_key, 0.957 + (0.001 if "post" in year else 0.0)
        )
        if hasattr(fj, "particleNet_TvsQCD"):
            top_score_pool    = ak.where(is_top_pool, fj.particleNet_TvsQCD, -1.0)
            leading_top_score = ak.max(top_score_pool, axis=1)
            any_top = leading_top_score > top_thresh
        else:
            any_top = ak.zeros_like(best_pt, dtype=bool)
        pass_notop = pass_leadPt & (~any_top)

        # nonHto4bFatJet for WZ veto:
        #   selFatJets with dR > 0.8 from Higgs candidate (standard)
        is_wz_pool = (
            base
            & (dR_from_best > 0.05)
            & (fj_pt > 250)
            & (fj_sd > 50)
            & (fj_sd < 220)
        )

        # WZ veto: leading WZvsQCD score in WZ pool > threshold
        wz_thresh  = {"2018": 0.9873, "2017": 0.9858}.get(self._trig_key, 0.9843)
        if hasattr(fj, "particleNet_WZvsQCD"):
            wz_score_pool    = ak.where(is_wz_pool, fj.particleNet_WZvsQCD, -1.0)
            leading_wz_score = ak.max(wz_score_pool, axis=1)
            any_wz = leading_wz_score > wz_thresh
        else:
            any_wz = ak.zeros_like(best_pt, dtype=bool)
        pass_nowz = pass_notop & (~any_wz)

        return (pass_candH, pass_leadPt, pass_notop, pass_nowz,
                best_pt, best_eta, best_phi, best_mass,
                best_x4b, best_massH_v2b, best_34massAa)

    # ------------------------------------------------------------------
    # Main process method
    # ------------------------------------------------------------------

    @staticmethod
    def _in_hem(eta, phi):
        """True for objects in the 2018 HEM15/16 affected region."""
        return (eta > -3.0) & (eta < -1.3) & (phi > -1.57) & (phi < -0.87)

    def process(self, events):
        year = self.year
        self.cnt["total"] = self.cnt.get("total", 0) + len(events)

        def _cf(label, n):
            self.cnt[label] = self.cnt.get(label, 0) + n
            print(f"  {label:<34s}: {n:>10,}")

        def _filt(mask, *arrs):
            """Filter events plus any number of parallel flat arrays."""
            return tuple(a[mask] for a in arrs)

        print(f"\n  {chr(9472)*48}")
        print(f"  Chunk total: {len(events):,}")

        # ── 1. run:ls ────────────────────────────────────────────────
        if self.lumi_mask is not None:
            events = events[self.lumi_mask(events.run, events.luminosityBlock)]
        _cf("run:ls", len(events))
        if len(events) == 0: return {}

        # ── 2. nPV ───────────────────────────────────────────────────
        events = events[events.PV.npvsGood > 0]
        _cf("nPV", len(events))
        if len(events) == 0: return {}

        # ── 3. METFilters ────────────────────────────────────────────
        try:
            events = events[self._noise_flags(events)]
        except Exception as exc:
            print(f"  WARNING noise flags: {exc}")
        _cf("METFilters", len(events))
        if len(events) == 0: return {}

        # ── 4+5+8+9. candH / leadPt / TopVeto / WZVeto ───────────────
        # All four masks come from _select_fatjet and are boolean arrays over
        # the SAME events array (cumulative: nowz⊆notop⊆leadPt⊆candH).
        # Apply them ONE AT A TIME so each cutflow step sees the right count.
        (pass_candH, pass_leadPt, pass_notop, pass_nowz,
         fj_pt, fj_eta, fj_phi, fj_mass,
         x4b, pnet_v2b, pnet_34) = self._select_fatjet(events)

        # 4. candH
        (events,
         fj_pt, fj_eta, fj_phi, fj_mass,
         x4b, pnet_v2b, pnet_34) = _filt(
            pass_candH,
            events, fj_pt, fj_eta, fj_phi, fj_mass,
            x4b, pnet_v2b, pnet_34)
        # Re-slice the remaining masks to match the now-smaller events array.
        # Since pass_leadPt[i] implies pass_candH[i], we extract the sub-mask:
        pass_leadPt_sub = pass_leadPt[pass_candH]
        pass_notop_sub  = pass_notop [pass_candH]
        pass_nowz_sub   = pass_nowz  [pass_candH]
        _cf("candH", len(events))
        if len(events) == 0: return {}

        # 5. leadingFatJetPt
        (events,
         fj_pt, fj_eta, fj_phi, fj_mass,
         x4b, pnet_v2b, pnet_34) = _filt(
            pass_leadPt_sub,
            events, fj_pt, fj_eta, fj_phi, fj_mass,
            x4b, pnet_v2b, pnet_34)
        pass_notop_sub = pass_notop_sub[pass_leadPt_sub]
        pass_nowz_sub  = pass_nowz_sub [pass_leadPt_sub]
        _cf("leadingFatJetPt", len(events))
        if len(events) == 0: return {}

        # ── 6. Trg_Combo ─────────────────────────────────────────────
        trig_mask, trig_prio = self._trigger_mask(events)
        (events, trig_prio,
         fj_pt, fj_eta, fj_phi, fj_mass,
         x4b, pnet_v2b, pnet_34) = _filt(
            trig_mask,
            events, trig_prio,
            fj_pt, fj_eta, fj_phi, fj_mass,
            x4b, pnet_v2b, pnet_34)
        pass_notop_sub = pass_notop_sub[trig_mask]
        pass_nowz_sub  = pass_nowz_sub [trig_mask]
        _cf("Trg_Combo_AK4AK8Jet_HT_VBF", len(events))
        if len(events) == 0: return {}

        # ── 7. nLeptonsTight ─────────────────────────────────────────
        lep_veto = self._muon_veto(events) & self._electron_veto(events)
        (events, trig_prio,
         fj_pt, fj_eta, fj_phi, fj_mass,
         x4b, pnet_v2b, pnet_34) = _filt(
            lep_veto,
            events, trig_prio,
            fj_pt, fj_eta, fj_phi, fj_mass,
            x4b, pnet_v2b, pnet_34)
        pass_notop_sub = pass_notop_sub[lep_veto]
        pass_nowz_sub  = pass_nowz_sub [lep_veto]
        _cf("nLeptonsTight", len(events))
        if len(events) == 0: return {}

        # ── 8. TopFatJetVeto ─────────────────────────────────────────
        (events, trig_prio,
         fj_pt, fj_eta, fj_phi, fj_mass,
         x4b, pnet_v2b, pnet_34) = _filt(
            pass_notop_sub,
            events, trig_prio,
            fj_pt, fj_eta, fj_phi, fj_mass,
            x4b, pnet_v2b, pnet_34)
        pass_nowz_sub = pass_nowz_sub[pass_notop_sub]
        _cf("TopFatJetVeto", len(events))
        if len(events) == 0: return {}

        # ── 9. VFatJetVjjVeto (WZ) ───────────────────────────────────
        (events, trig_prio,
         fj_pt, fj_eta, fj_phi, fj_mass,
         x4b, pnet_v2b, pnet_34) = _filt(
            pass_nowz_sub,
            events, trig_prio,
            fj_pt, fj_eta, fj_phi, fj_mass,
            x4b, pnet_v2b, pnet_34)
        _cf("VFatJetVjjVeto", len(events))
        if len(events) == 0: return {}

        # ── 10. MetZvvVeto ────────────────────────────────────────────
        # Reference uses MET_T1_pt (Type-1 corrected), not raw MET_pt
        met_pt = events.MET.T1_pt if hasattr(events.MET, "T1_pt") else events.MET.pt
        met_ok = met_pt <= 200
        (events, trig_prio,
         fj_pt, fj_eta, fj_phi, fj_mass,
         x4b, pnet_v2b, pnet_34) = _filt(
            met_ok,
            events, trig_prio,
            fj_pt, fj_eta, fj_phi, fj_mass,
            x4b, pnet_v2b, pnet_34)
        _cf("MetZvvVeto", len(events))
        if len(events) == 0: return {}

        # ── 11. AK4 jets ─────────────────────────────────────────────
        # Mirrors reference selectAK4Jets + b-veto logic:
        #
        # b-veto keeps a jet if:
        #   btagDeepFlavB <= WP  (not b-tagged)
        #   OR |eta| >= 2.4      (forward jets always kept regardless of b-tag)
        #
        # Then VBF selection uses leading 2 of these b-vetoed jets (no separate
        # VBF veto step — just one selection on b-veto jets).
        j      = events.Jet
        j_pt   = self._pt_field(j)
        j_mass = self._mass_field(j)
        bWP    = self.btagWPM

        j_ok = (
            (j_pt > 30)
            & (j.jetId == 6)
            & ((j_pt > 50) | (j.puId >= 4))
            & (self._delta_r(j.eta, j.phi, fj_eta[:, None], fj_phi[:, None]) > 0.8)
        )
        jets_sel = j    [j_ok]
        jsel_pt  = j_pt [j_ok]
        jsel_m   = j_mass[j_ok]

        # b-jet = central (|eta|<2.4) AND b-tagged; forward jets always pass veto
        is_bjet  = (abs(jets_sel.eta) < 2.4) & (jets_sel.btagDeepFlavB > bWP)
        has_bjet = ak.any(is_bjet, axis=1)

        # BJetVeto: apply as hard filter (matching reference cutflow order)
        _cf("BJetVeto", int(ak.sum(~has_bjet)))
        (events, trig_prio,
         fj_pt, fj_eta, fj_phi, fj_mass,
         x4b, pnet_v2b, pnet_34,
         jets_sel, jsel_pt, jsel_m, is_bjet) = _filt(
            ~has_bjet,
            events, trig_prio,
            fj_pt, fj_eta, fj_phi, fj_mass,
            x4b, pnet_v2b, pnet_34,
            jets_sel, jsel_pt, jsel_m, is_bjet)
        if len(events) == 0: return {}

        # b-vetoed jets: NOT b-tagged OR forward (matches Sel_Ak4Jets_bveto)
        bveto_jets  = jets_sel[~is_bjet]
        bveto_pt    = jsel_pt [~is_bjet]
        bveto_mass  = jsel_m  [~is_bjet]

        def _lead2(arr, pad):
            return ak.fill_none(ak.pad_none(arr, 2, clip=True), pad)

        # Apply n_bveto >= 2 before computing kinematics (matching reference
        # which uses ak.mask to set events with < 2 jets to None)
        n_bveto = ak.num(bveto_pt)
        (events, trig_prio,
         fj_pt, fj_eta, fj_phi, fj_mass,
         x4b, pnet_v2b, pnet_34,
         bveto_jets, bveto_pt, bveto_mass) = _filt(
            n_bveto >= 2,
            events, trig_prio,
            fj_pt, fj_eta, fj_phi, fj_mass,
            x4b, pnet_v2b, pnet_34,
            bveto_jets, bveto_pt, bveto_mass)
        if len(events) == 0: return {}

        # VBF selection: leading 2 of b-vetoed jets (pt-sorted)
        '''
        sort_bv  = ak.argsort(bveto_pt, ascending=False)
        q1_pt    = ak.flatten(_lead2(bveto_pt       [sort_bv], 0.0)[:, 0:1], axis=1)
        q1_eta   = ak.flatten(_lead2(bveto_jets.eta [sort_bv], 0.0)[:, 0:1], axis=1)
        q1_phi   = ak.flatten(_lead2(bveto_jets.phi [sort_bv], 0.0)[:, 0:1], axis=1)
        q1_mass  = ak.flatten(_lead2(bveto_mass     [sort_bv], 0.0)[:, 0:1], axis=1)
        q2_pt    = ak.flatten(_lead2(bveto_pt       [sort_bv], 0.0)[:, 1:2], axis=1)
        q2_eta   = ak.flatten(_lead2(bveto_jets.eta [sort_bv], 0.0)[:, 1:2], axis=1)
        q2_phi   = ak.flatten(_lead2(bveto_jets.phi [sort_bv], 0.0)[:, 1:2], axis=1)
        q2_mass  = ak.flatten(_lead2(bveto_mass     [sort_bv], 0.0)[:, 1:2], axis=1)
        deta_sel = abs(q1_eta - q2_eta)
        mjj_sel  = self._dijet_mass(q1_pt, q1_eta, q1_phi, q1_mass,
                                    q2_pt, q2_eta, q2_phi, q2_mass)
        '''
        # VBF selection: leading 2 of b-vetoed jets (as stored in file, no sorting)
        q1_pt    = ak.flatten(_lead2(bveto_pt, 0.0)[:, 0:1], axis=1)
        q1_eta   = ak.flatten(_lead2(bveto_jets.eta, 0.0)[:, 0:1], axis=1)
        q1_phi   = ak.flatten(_lead2(bveto_jets.phi, 0.0)[:, 0:1], axis=1)
        q1_mass  = ak.flatten(_lead2(bveto_mass, 0.0)[:, 0:1], axis=1)
        q2_pt    = ak.flatten(_lead2(bveto_pt, 0.0)[:, 1:2], axis=1)
        q2_eta   = ak.flatten(_lead2(bveto_jets.eta, 0.0)[:, 1:2], axis=1)
        q2_phi   = ak.flatten(_lead2(bveto_jets.phi, 0.0)[:, 1:2], axis=1)
        q2_mass  = ak.flatten(_lead2(bveto_mass, 0.0)[:, 1:2], axis=1)
        deta_sel = abs(q1_eta - q2_eta)
        mjj_sel  = self._dijet_mass(q1_pt, q1_eta, q1_phi, q1_mass,
                                    q2_pt, q2_eta, q2_phi, q2_mass)
        # nak4jets = dEta > 2.2 AND mjj > 450 (BJetVeto already applied above)
        pass_nak4 = (deta_sel > 2.2) & (mjj_sel > 450)
        (events, trig_prio,
         fj_pt, fj_eta, fj_phi, fj_mass,
         x4b, pnet_v2b, pnet_34,
         deta_sel, mjj_sel,
         q1_pt, q1_eta, q1_phi, q1_mass,
         q2_pt, q2_eta, q2_phi, q2_mass) = _filt(
            pass_nak4,
            events, trig_prio,
            fj_pt, fj_eta, fj_phi, fj_mass,
            x4b, pnet_v2b, pnet_34,
            deta_sel, mjj_sel,
            q1_pt, q1_eta, q1_phi, q1_mass,
            q2_pt, q2_eta, q2_phi, q2_mass)
        _cf("nak4jets", len(events))
        if len(events) == 0: return {}

        # ── 12. 2018HEM1516Issue ──────────────────────────────────────
        # Only applies to runs affected by HEM15/16 issue: [319077, 325175]
        # Different eta/phi boundaries for fat jet vs AK4 jets:
        #   FatJet: eta < -1.1 AND |phi + 1.22| < 0.55
        #   AK4:    -3.2 < eta < -1.2 AND |phi + 1.22| < 0.45
        if "2018" in year:
            is_affected_run = (
                (events.run >= 319077) & (events.run <= 325175)
            )
            hem_fj = (fj_eta < -1.1) & (abs(fj_phi + 1.22) < 0.55)
            hem_q1 = (q1_eta > -3.2) & (q1_eta < -1.2) & (abs(q1_phi + 1.22) < 0.45)
            hem_q2 = (q2_eta > -3.2) & (q2_eta < -1.2) & (abs(q2_phi + 1.22) < 0.45)
            hem_ok = ~(is_affected_run & (hem_fj | hem_q1 | hem_q2))
            (events, trig_prio,
             fj_pt, fj_eta, fj_phi, fj_mass,
             x4b, pnet_v2b, pnet_34,
             deta_sel, mjj_sel,
             q1_pt, q1_eta, q1_phi, q1_mass,
             q2_pt, q2_eta, q2_phi, q2_mass) = _filt(
                hem_ok,
                events, trig_prio,
                fj_pt, fj_eta, fj_phi, fj_mass,
                x4b, pnet_v2b, pnet_34,
                deta_sel, mjj_sel,
                q1_pt, q1_eta, q1_phi, q1_mass,
                q2_pt, q2_eta, q2_phi, q2_mass)
            _cf("2018HEM1516Issue", len(events))
            if len(events) == 0: return {}

        # ── 13. VBF sub-category ─────────────────────────────────────
        '''
        is_hi    = (deta_sel > 3.0) & (mjj_sel > 900)
        is_hi_pt = fj_pt >= 400
        vbf_cat_idx = ak.where(
            is_hi & is_hi_pt,        np.int32(3),
            ak.where(is_hi & (~is_hi_pt), np.int32(2),
            ak.where((~is_hi) & is_hi_pt, np.int32(1), np.int32(0))))
        '''
        # ── 13. VBF sub-category ─────────────────────────────────────
        is_hi    = (deta_sel > 3.0) & (mjj_sel > 900)
        is_lo    = (deta_sel <= 3.0) | (mjj_sel <= 900)  # This matches the new reference logic
        is_lo_pt = fj_pt <= 400
        is_hi_pt = fj_pt > 400
        vbf_cat_idx = ak.where(
            is_hi & is_hi_pt,        np.int32(3),   # VBFjjHiPtHi (pt > 400)
            ak.where(is_hi & is_lo_pt, np.int32(2),  # VBFjjHiPtLo (pt <= 400)
                     ak.where(is_lo & is_hi_pt, np.int32(1),  # VBFjjLoPtHi (pt > 400)
                              ak.where(is_lo & is_lo_pt, np.int32(0),  # VBFjjLoPtLo (pt <= 400)
                                       np.int32(-1)))))
            
        for ci, cn in enumerate(VBF_CATEGORIES):
            _cf(f"  {cn}", int(ak.sum(vbf_cat_idx == ci)))

        # ── 14. X4b regions ──────────────────────────────────────────
        is_SR      = x4b > 0.96
        is_SB      = x4b > 0.84
        is_SB_only = (x4b > 0.84) & (x4b <= 0.96)
        _cf("X4bSR  (score>0.96)", int(ak.sum(is_SR)))
        _cf("X4bSB  (score>0.84)", int(ak.sum(is_SB)))

        # ── 15. Fill output buffers ───────────────────────────────────
        to_np = lambda a: np.asarray(ak.to_numpy(a))
        run_np  = to_np(events.run);          ls_np  = to_np(events.luminosityBlock)
        evt_np  = to_np(events.event)
        fpt_np  = to_np(fj_pt);               feta_np = to_np(fj_eta)
        fphi_np = to_np(fj_phi);              fmas_np = to_np(fj_mass)
        q1pt_np = to_np(q1_pt);              q1et_np = to_np(q1_eta)
        q1ph_np = to_np(q1_phi);             q1ms_np = to_np(q1_mass)
        q2pt_np = to_np(q2_pt);              q2et_np = to_np(q2_eta)
        q2ph_np = to_np(q2_phi);             q2ms_np = to_np(q2_mass)
        det_np  = to_np(deta_sel);            mjj_np  = to_np(mjj_sel)
        x4b_np  = to_np(x4b);                pv2b_np = to_np(pnet_v2b)
        p34_np  = to_np(pnet_34)
        prio_np = to_np(trig_prio).astype(np.int32)
        tf_np   = np.ones(len(events), dtype=np.int32)
        cat_np  = to_np(vbf_cat_idx)
        sr_np   = to_np(is_SR).astype(bool);  sb_np  = to_np(is_SB).astype(bool)
        sbo_np  = to_np(is_SB_only).astype(bool)

        for region_name, rm in [("X4bSR", sr_np), ("X4bSB", sb_np), ("X4bSB_only", sbo_np)]:
            for cat_idx, cat_name in enumerate(VBF_CATEGORIES):
                sel = rm & (cat_np == cat_idx)
                if not np.any(sel): continue
                buf = self._out[f"{cat_name}_{region_name}"]
                buf["run"]             .extend(run_np [sel].tolist())
                buf["luminosityBlock"] .extend(ls_np  [sel].tolist())
                buf["event"]           .extend(evt_np [sel].tolist())
                buf["fatjet_pt"]       .extend(fpt_np [sel].tolist())
                buf["fatjet_eta"]      .extend(feta_np[sel].tolist())
                buf["fatjet_phi"]      .extend(fphi_np[sel].tolist())
                buf["fatjet_mass"]     .extend(fmas_np[sel].tolist())
                buf["vbfjet1_pt"]      .extend(q1pt_np[sel].tolist())
                buf["vbfjet1_eta"]     .extend(q1et_np[sel].tolist())
                buf["vbfjet1_phi"]     .extend(q1ph_np[sel].tolist())
                buf["vbfjet1_mass"]    .extend(q1ms_np[sel].tolist())
                buf["vbfjet2_pt"]      .extend(q2pt_np[sel].tolist())
                buf["vbfjet2_eta"]     .extend(q2et_np[sel].tolist())
                buf["vbfjet2_phi"]     .extend(q2ph_np[sel].tolist())
                buf["vbfjet2_mass"]    .extend(q2ms_np[sel].tolist())
                buf["vbfjj_deta"]      .extend(det_np [sel].tolist())
                buf["vbfjj_mass"]      .extend(mjj_np [sel].tolist())
                buf["x4b_score"]       .extend(x4b_np [sel].tolist())
                buf["pnet_massH_v2b"]  .extend(pv2b_np[sel].tolist())
                buf["pnet_34massAa"]   .extend(p34_np [sel].tolist())
                buf["trigger_fired"]   .extend(tf_np  [sel].tolist())
                buf["trigger_priority"].extend(prio_np[sel].tolist())
        return {}

    def postprocess(self, accumulator):
        return accumulator

    # ------------------------------------------------------------------
    # Write output ROOT file (call once after all files are processed)
    # ------------------------------------------------------------------

    def flush_output(self):
        with uproot.recreate(self.output_file) as f:
            for tree_key, branches in self._out.items():
                if not branches["run"]:
                    continue   # skip empty trees
                data = {k: np.array(v) for k, v in branches.items()}
                f[tree_key] = data
        print(f"\nOutput written to: {self.output_file}")

    def print_summary(self, dataset_name):
        CUTFLOW_ORDER = [
            "total", "run:ls", "nPV", "METFilters",
            "candH", "leadingFatJetPt",
            "Trg_Combo_AK4AK8Jet_HT_VBF",
            "nLeptonsTight", "TopFatJetVeto", "VFatJetVjjVeto",
            "MetZvvVeto", "BJetVeto", "nak4jets", "2018HEM1516Issue",
        ] + [f"  {c}" for c in VBF_CATEGORIES] + [
            "X4bSR  (score>0.96)", "X4bSB  (score>0.84)",
        ]
        print(f"\n{'='*70}")
        print(f"Cutflow: {dataset_name}")
        print(f"{'='*70}")
        print(f"  {'Cut':<34s}  {'Events':>12s}")
        print(f"  {'-'*34}  {'-'*12}")
        for k in CUTFLOW_ORDER:
            if k in self.cnt:
                print(f"  {k:<34s}  {self.cnt[k]:>12,}")
        print(f"\n{'='*70}")
        print("Tree counts:")
        for cat in VBF_CATEGORIES:
            print(f"\n  {cat}:")
            for sel in X4B_SELECTIONS:
                n = len(self._out[f"{cat}_{sel}"]["run"])
                print(f"    {sel}: {n}")



# ──────────────────────────────────────────────────────────────────────────────
# Per-file worker (runs in a subprocess when using multiprocessing)
# ──────────────────────────────────────────────────────────────────────────────

def _process_one_file(args):
    """
    Process a single ROOT file and return the output buffers + counters.
    This function is called in a worker process.
    """
    fpath, year, primary_dataset, lumi_mask_path, chunk_size, verbose = args

    # Re-create LumiMask inside the worker (not picklable)
    lumi_mask = LumiMask(lumi_mask_path) if lumi_mask_path else None

    proc = VBFx4bProcessor(year, primary_dataset, lumi_mask, output_file=None)

    fname = os.path.basename(fpath)
    t0 = time.time()

    try:
        # Use uproot.iterate for memory-efficient chunked reading
        for batch in uproot.iterate(
            {fpath: "Events"},
            filter_name=None,          # read all branches (NanoEvents schema needs them)
            step_size=chunk_size,
            library="ak",
        ):
            # Wrap the raw awkward batch in NanoEvents so field access works
            events = NanoEventsFactory.from_root(
                fpath,
                schemaclass=NanoAODSchema,
                entry_start=int(batch["run"].layout.length) if hasattr(batch["run"].layout, "length") else 0,
            ).events()
            # NOTE: uproot.iterate gives raw ak arrays; NanoEventsFactory gives
            # proper cross-referenced objects. We use NanoEventsFactory per-chunk
            # via entry_start/entry_stop below instead.
            break  # escape the test loop

    except Exception:
        pass  # fall through to the correct approach below

    # Correct approach: NanoEventsFactory with entry_start/entry_stop chunking
    try:
        f_root = uproot.open(fpath)
        tree   = f_root["Events"]
        n_total = tree.num_entries
    except Exception as exc:
        return None, f"Cannot open {fpath}: {exc}"

    n_chunks = max(1, (n_total + chunk_size - 1) // chunk_size)
    desc = fname[:40]

    if verbose and HAS_TQDM:
        pbar = tqdm(total=n_total, desc=desc, unit="ev", leave=False)
    else:
        pbar = None

    for i_chunk in range(n_chunks):
        entry_start = i_chunk * chunk_size
        entry_stop  = min(entry_start + chunk_size, n_total)
        try:
            events = NanoEventsFactory.from_root(
                fpath,
                schemaclass=NanoAODSchema,
                entry_start=entry_start,
                entry_stop=entry_stop,
            ).events()
            proc.process(events)
        except Exception as exc:
            print(f"  WARNING chunk {i_chunk} of {fname}: {exc}")
        if pbar:
            pbar.update(entry_stop - entry_start)

    if pbar:
        pbar.close()

    elapsed = time.time() - t0
    rate = n_total / elapsed if elapsed > 0 else 0
    if verbose:
        print(f"  {fname}: {n_total} events in {elapsed:.1f}s  ({rate/1e3:.1f}k ev/s)")

    return proc._out, proc.cnt


def _merge_outputs(all_outs, all_cnts):
    """Merge output buffers and counters from multiple workers."""
    merged_out = {
        f"{c}_{s}": {k: [] for k in (
            "run", "luminosityBlock", "event",
            "fatjet_pt", "fatjet_eta", "fatjet_phi", "fatjet_mass",
            "vbfjet1_pt", "vbfjet1_eta", "vbfjet1_phi", "vbfjet1_mass",
            "vbfjet2_pt", "vbfjet2_eta", "vbfjet2_phi", "vbfjet2_mass",
            "vbfjj_deta", "vbfjj_mass", "x4b_score",
            "pnet_massH_v2b", "pnet_34massAa",
            "trigger_fired", "trigger_priority",
        )}
        for c in VBF_CATEGORIES for s in X4B_SELECTIONS
    }
    merged_cnt = {k: 0 for k in all_cnts[0]}

    for out, cnt in zip(all_outs, all_cnts):
        for key in merged_out:
            for branch in merged_out[key]:
                merged_out[key][branch].extend(out[key][branch])
        for k in merged_cnt:
            merged_cnt[k] += cnt[k]

    return merged_out, merged_cnt


# ──────────────────────────────────────────────────────────────────────────────
# Dataset processing
# ──────────────────────────────────────────────────────────────────────────────

def process_dataset(dataset_name, output_name=None, n_workers=4, chunk_size=200_000, verbose=True):
    year = extract_year_from_dataset(dataset_name)
    print(f"\n{'='*70}")
    print(f"Processing : {dataset_name}")
    print(f"Year       : {year}")
    print(f"Workers    : {n_workers}   Chunk size: {chunk_size:,}")

    # Golden JSON — write to temp file for pickling across workers
    golden_json_str = download_golden_json(year)
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False)
    tmp.write(golden_json_str)
    tmp.close()
    lumi_mask_path = tmp.name

    # Input files
    input_files = get_file_list(dataset_name, year)
    if not input_files:
        os.unlink(lumi_mask_path)
        return False
    input_files = [f for f in input_files if os.path.exists(f)]
    if not input_files:
        print("Error: no files found on disk")
        os.unlink(lumi_mask_path)
        return False

    primary_dataset = extract_primary_dataset(input_files)
    print(f"Primary dataset: {primary_dataset}")
    print(f"Found {len(input_files)} input file(s)")

    if output_name is None:
        output_name = f"{dataset_name}_output.root"

    # Build args for each file
    file_args = [
        (fpath, year, primary_dataset, lumi_mask_path, chunk_size, verbose)
        for fpath in input_files
    ]

    t_start = time.time()
    all_outs = []
    all_cnts = []

    actual_workers = min(n_workers, len(input_files))

    if actual_workers <= 1:
        # Single-process path (easier to debug)
        for fa in file_args:
            out, cnt = _process_one_file(fa)
            if out is not None:
                all_outs.append(out)
                all_cnts.append(cnt)
            else:
                print(f"  WARNING: {cnt}")
    else:
        # Multi-process path
        with ProcessPoolExecutor(max_workers=actual_workers) as pool:
            futures = {pool.submit(_process_one_file, fa): fa[0] for fa in file_args}
            if HAS_TQDM:
                pbar = tqdm(total=len(futures), desc="Files", unit="file")
            for fut in as_completed(futures):
                fpath = futures[fut]
                try:
                    out, cnt = fut.result()
                    if out is not None:
                        all_outs.append(out)
                        all_cnts.append(cnt)
                    else:
                        print(f"  WARNING {os.path.basename(fpath)}: {cnt}")
                except Exception as exc:
                    print(f"  ERROR {os.path.basename(fpath)}: {exc}")
                if HAS_TQDM:
                    pbar.update(1)
            if HAS_TQDM:
                pbar.close()

    os.unlink(lumi_mask_path)

    if not all_outs:
        print("No output produced.")
        return False

    # Merge all worker outputs
    merged_out, merged_cnt = _merge_outputs(all_outs, all_cnts)

    # Write ROOT file
    with uproot.recreate(output_name) as f:
        for tree_key, branches in merged_out.items():
            if not branches["run"]:
                continue
            f[tree_key] = {k: np.array(v) for k, v in branches.items()}

    elapsed = time.time() - t_start

    # Summary
    print(f"\n{'='*70}")
    print(f"Summary for {dataset_name}:")
    print(f"{'='*70}")
    CUTFLOW_ORDER = [
        "total", "run:ls", "nPV", "METFilters",
        "candH", "leadingFatJetPt", "Trg_Combo_AK4AK8Jet_HT_VBF",
        "nLeptonsTight", "TopFatJetVeto", "VFatJetVjjVeto",
        "MetZvvVeto", "BJetVeto", "nak4jets", "2018HEM1516Issue",
    ] + [f"  {c}" for c in VBF_CATEGORIES] + ["X4bSR  (score>0.96)", "X4bSB  (score>0.84)"]
    print(f"  {'Cut':<34s}  {'Events':>12s}")
    print(f"  {'-'*34}  {'-'*12}")
    for k in CUTFLOW_ORDER:
        if k in merged_cnt:
            print(f"  {k:<34s}  {merged_cnt[k]:>12,}")
    print(f"Wall time              : {elapsed:.1f}s  ({merged_cnt['total']/elapsed/1e3:.1f}k ev/s)")
    print("\nTree counts:")
    for cat in VBF_CATEGORIES:
        print(f"\n  {cat}:")
        for sel in X4B_SELECTIONS:
            key = f"{cat}_{sel}"
            n = len(merged_out[key]["run"])
            print(f"    {sel}: {n}")
    print(f"\nOutput written to: {output_name}")
    return True


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Process VBF events with X4b selections (Coffea version)"
    )
    parser.add_argument("datasets", nargs="+", help="Dataset name(s) to process")
    parser.add_argument("--output",     "-o", default=None,
                        help="Output ROOT file name (default: {dataset}_output.root)")
    parser.add_argument("--workers",    "-j", type=int, default=4,
                        help="Number of parallel worker processes (default: 4)")
    parser.add_argument("--chunk-size", "-c", type=int, default=200_000,
                        help="Events per chunk per file (default: 200000)")
    parser.add_argument("--quiet",      "-q", action="store_true",
                        help="Suppress per-file progress output")
    args = parser.parse_args()

    print(f"\n{'='*70}")
    print(f"Processing {len(args.datasets)} dataset(s)")
    print(f"{'='*70}")

    for ds in args.datasets:
        process_dataset(
            ds,
            output_name=args.output,
            n_workers=args.workers,
            chunk_size=args.chunk_size,
            verbose=not args.quiet,
        )

    print(f"\n{'='*70}")
    print("All done!")
    print(f"{'='*70}\n")


if __name__ == "__main__":
    main()
