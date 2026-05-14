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
            "AK8PFJet500", "PFHT1050", "PFJet500",
            "AK8PFHT800_TrimMass50", "AK8PFJet400_TrimMass30",
            "AK8PFJet360_TrimMass30",
        ],
        "BTagCSV": [
            "PFHT380_SixPFJet32_DoublePFBTagCSV_2p2",
            "PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2",
            "PFHT430_SixPFJet40_PFBTagCSV_1p5",
            "AK8PFHT750_TrimMass50", "AK8PFJet380_TrimMass30",
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
        ],
        "BTagCSV": [
            "AK8DiPFJet250_200_TrimMass30_BTagCSV_p20",
            "AK8DiPFJet280_200_TrimMass30_BTagCSV_p20",
            "PFHT400_SixJet30_DoubleBTagCSV_p056",
            "PFHT450_SixJet40_BTagCSV_p056",
            "AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20",
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
        self.btagWPM = {"2018": 0.2783, "2017": 0.3040}.get(
            self._trig_key, 0.2489
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

        # Detailed cutflow counters
        self.cnt = {k: 0 for k in (
            "total", "pass_golden", "pass_pv", "pass_noise",
            "pass_trigger", "pass_fatH", "pass_lepton",
            "pass_met", "pass_VBF",
        )}

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
            mask = fired_jht | fired_btag
            prio = ak.where(fired_jht, np.int32(1), np.int32(2))
            prio = ak.where(mask, prio, np.int32(0))
        elif pd == "BTagCSV":
            mask = fired_btag & (~fired_jht)
            prio = ak.where(mask, np.int32(2), np.int32(0))
        else:
            mask = fired_jht | fired_btag
            prio = ak.zeros_like(events.run, dtype=np.int32)

        return mask, prio

    # ------------------------------------------------------------------
    # Selection helpers
    # ------------------------------------------------------------------

    def _noise_flags(self, events):
        F = events.Flag
        base = (
            F.goodVertices
            & F.globalSuperTightHalo2016Filter
            & F.HBHENoiseFilter
            & F.HBHENoiseIsoFilter
            & F.eeBadScFilter
            & F.BadPFMuonFilter
            & F.BadPFMuonDzFilter
            & F.EcalDeadCellTriggerPrimitiveFilter
        )
        if "2016" in self.year:
            base = base & F.ecalBadCalibFilter
        return base

    def _muon_veto(self, events):
        mu = events.Muon
        pt_thr = 26.0 + (3.0 if "2017" in self.year else 0.0)
        good = (
            (mu.pt > pt_thr)
            & (abs(mu.eta) < 2.4)
            & (abs(mu.dxy) < 0.02)
            & (abs(mu.dz)  < 0.10)
            & (mu.miniPFRelIso_all < 0.10)
            & (mu.mediumId | ((mu.pt > 200) & (mu.highPtId > 0)))
        )
        return ~ak.any(good, axis=1)

    def _electron_veto(self, events):
        ele = events.Electron
        pt_thr = 30.0 + (5.0 if "2016" not in self.year else 0.0)
        in_crack = (abs(ele.eta) > 1.44) & (abs(ele.eta) < 1.57)
        wp90 = ele.mvaFall17V2Iso_WP90
        wp80 = ele.mvaFall17V2Iso_WP80
        heep = ele.cutBased_HEEP
        good = (
            (ele.pt > pt_thr)
            & (abs(ele.eta) < 2.5)
            & (~in_crack)
            & wp90
            & (wp80 | ((ele.pt > 35) & (heep > 0)))
            & (abs(ele.dxy) < 0.02)
            & (abs(ele.dz)  < 0.10)
        )
        return ~ak.any(good, axis=1)

    def _select_fatjet(self, events):
        """
        Select the best Higgs-candidate fat jet per event.
        Returns (pass_mask, pt, eta, phi, mass, x4b, pnet_v2b, pnet_34) as flat arrays.
        All returned arrays are per-event (length = len(events)).
        Events where no candidate is found have pass_mask=False; their kinematic
        values are undefined (but will be masked out before use).
        """
        fj   = events.FatJet
        year = self.year

        fj_pt   = self._pt_field(fj)
        fj_mass = self._mass_field(fj)
        fj_sd   = self._msoftdrop_field(fj)

        # Basic kinematic / quality cuts
        base = (
            (fj_pt > 250)
            & (abs(fj.eta) < 2.4)
            & (fj.jetId == 6)
            & (fj_sd > 20)
        )

        # Xbb tagger
        xbb_ok = fj.particleNetMD_XbbvsQCD > 0.75

        # X4b score
        if hasattr(fj, "PNet_X4b_v2a_Haa4b_score") and hasattr(fj, "PNet_X4b_v2b_Haa4b_score"):
            x4b_scores = 0.5 * (fj.PNet_X4b_v2a_Haa4b_score + fj.PNet_X4b_v2b_Haa4b_score)
        elif hasattr(fj, "X4b_score"):
            x4b_scores = fj.X4b_score
        else:
            x4b_scores = ak.full_like(fj_pt, -999.0)

        # Candidate mask: must pass base + Xbb
        cand_mask = base & xbb_ok

        # Score for argmax: -999 for non-candidates
        score_for_max = ak.where(cand_mask, x4b_scores, -999.0)

        # Top tagger threshold
        top_thresh = {"2018": 0.970, "2017": 0.970}.get(
            self._trig_key, 0.957 + (0.001 if "post" in year else 0.0)
        )
        if hasattr(fj, "particleNet_TvsQCD"):
            top_tagged = (fj.particleNet_TvsQCD > top_thresh) & (fj_pt > 300)
        else:
            top_tagged = ak.zeros_like(fj_pt, dtype=bool)

        wz_thresh = {"2018": 0.9873, "2017": 0.9858}.get(self._trig_key, 0.9843)
        if hasattr(fj, "particleNet_WZvsQCD"):
            wz_tagged = fj.particleNet_WZvsQCD > wz_thresh
        else:
            wz_tagged = ak.zeros_like(fj_pt, dtype=bool)

        # ── Best candidate per event via argmax ──────────────────────
        # argmax returns an OptionType when keepdims=False on empty arrays;
        # use keepdims=True then index carefully.
        has_any_cand = ak.any(cand_mask, axis=1)

        # For events with no candidates, pad with a dummy 0 to avoid errors;
        # has_any_cand will mask those out later.
        n_fj = ak.num(fj_pt)
        safe_score = ak.where(n_fj > 0, score_for_max,
                              ak.Array([[0.0]] * len(events)))  # never reached
        best_idx_flat = ak.argmax(safe_score, axis=1)           # shape (N,) int

        # Gather per-event best-jet quantities using advanced indexing
        def _gather(arr):
            """Pick arr[i][best_idx_flat[i]] for each event i."""
            return arr[ak.local_index(arr) == best_idx_flat[:, None]]
            # Result is a ragged array with exactly 1 entry per event where
            # a match exists; flatten to get a flat array.

        best_x4b   = ak.flatten(_gather(score_for_max),  axis=None)
        best_pt    = ak.flatten(_gather(fj_pt),           axis=None)
        best_eta   = ak.flatten(_gather(fj.eta),          axis=None)
        best_phi   = ak.flatten(_gather(fj.phi),          axis=None)
        best_mass  = ak.flatten(_gather(fj_mass),         axis=None)

        if hasattr(fj, "PNet_massH_v2b"):
            best_massH_v2b = ak.flatten(_gather(fj.PNet_massH_v2b), axis=None)
        else:
            best_massH_v2b = ak.full_like(best_pt, -999.0)

        if hasattr(fj, "PNet_34massAa"):
            best_34massAa = ak.flatten(_gather(fj.PNet_34massAa), axis=None)
        else:
            best_34massAa = ak.full_like(best_pt, -999.0)

        # Top/WZ veto: any jet OTHER than the best candidate
        is_best_jet = ak.local_index(fj_pt) == best_idx_flat[:, None]
        any_top = ak.any(top_tagged & (~is_best_jet), axis=1)
        any_wz  = ak.any(wz_tagged  & (~is_best_jet), axis=1)
        veto    = any_top | any_wz

        pass_fatH = has_any_cand & (best_x4b > -999.0) & (~veto)

        return (pass_fatH,
                best_pt, best_eta, best_phi, best_mass,
                best_x4b, best_massH_v2b, best_34massAa)

    # ------------------------------------------------------------------
    # Main process method
    # ------------------------------------------------------------------

    def process(self, events):
        year = self.year
        n0 = len(events)
        self.cnt["total"] += n0

        # ── 1. Golden JSON ────────────────────────────────────────────
        if self.lumi_mask is not None:
            lumi_ok = self.lumi_mask(events.run, events.luminosityBlock)
            events = events[lumi_ok]
        self.cnt["pass_golden"] += len(events)
        if len(events) == 0:
            return {}

        # ── 2. Good PV ───────────────────────────────────────────────
        events = events[events.PV.npvsGood > 0]
        self.cnt["pass_pv"] += len(events)
        if len(events) == 0:
            return {}

        # ── 3. Noise flags ───────────────────────────────────────────
        try:
            events = events[self._noise_flags(events)]
        except Exception as exc:
            pass  # non-fatal; some flags may be absent in certain eras
        self.cnt["pass_noise"] += len(events)
        if len(events) == 0:
            return {}

        # ── 4. Trigger ───────────────────────────────────────────────
        trig_mask, trig_prio = self._trigger_mask(events)
        events    = events   [trig_mask]
        trig_prio = trig_prio[trig_mask]
        self.cnt["pass_trigger"] += len(events)
        if len(events) == 0:
            return {}

        # ── 5. Fat-jet selection ─────────────────────────────────────
        (pass_fatH,
         fj_pt, fj_eta, fj_phi, fj_mass,
         x4b, pnet_v2b, pnet_34) = self._select_fatjet(events)

        events    = events   [pass_fatH]
        trig_prio = trig_prio[pass_fatH]
        fj_pt     = fj_pt    [pass_fatH]
        fj_eta    = fj_eta   [pass_fatH]
        fj_phi    = fj_phi   [pass_fatH]
        fj_mass   = fj_mass  [pass_fatH]
        x4b       = x4b      [pass_fatH]
        pnet_v2b  = pnet_v2b [pass_fatH]
        pnet_34   = pnet_34  [pass_fatH]
        self.cnt["pass_fatH"] += len(events)
        if len(events) == 0:
            return {}

        # ── 6. Lepton vetoes ─────────────────────────────────────────
        mu_veto  = self._muon_veto    (events)
        ele_veto = self._electron_veto(events)
        lep_veto = mu_veto & ele_veto

        events    = events   [lep_veto]
        trig_prio = trig_prio[lep_veto]
        fj_pt     = fj_pt    [lep_veto]
        fj_eta    = fj_eta   [lep_veto]
        fj_phi    = fj_phi   [lep_veto]
        fj_mass   = fj_mass  [lep_veto]
        x4b       = x4b      [lep_veto]
        pnet_v2b  = pnet_v2b [lep_veto]
        pnet_34   = pnet_34  [lep_veto]
        self.cnt["pass_lepton"] += len(events)
        if len(events) == 0:
            return {}

        # ── 7. MET veto ──────────────────────────────────────────────
        met_ok = events.MET.pt <= 200
        events    = events   [met_ok]
        trig_prio = trig_prio[met_ok]
        fj_pt     = fj_pt    [met_ok]
        fj_eta    = fj_eta   [met_ok]
        fj_phi    = fj_phi   [met_ok]
        fj_mass   = fj_mass  [met_ok]
        x4b       = x4b      [met_ok]
        pnet_v2b  = pnet_v2b [met_ok]
        pnet_34   = pnet_34  [met_ok]
        self.cnt["pass_met"] += len(events)
        if len(events) == 0:
            return {}

        # ── 8. AK4 jet selection ─────────────────────────────────────
        j      = events.Jet
        j_pt   = self._pt_field(j)
        j_mass = self._mass_field(j)
        bWP    = self.btagWPM

        # Basic cuts
        j_ok = (
            (j_pt > 30)
            & (j.jetId == 6)
            & ((j_pt > 50) | (j.puId >= 4))
        )

        # dR > 0.8 from fat-jet candidate
        # fj_eta/phi are flat (N,); j.eta is ragged (N, njet)
        # Use awkward broadcasting: fj_eta[:,None] broadcasts over jet axis
        dR_jfj = self._delta_r(j.eta, j.phi,
                                fj_eta[:, None], fj_phi[:, None])
        j_ok = j_ok & (dR_jfj > 0.8)

        jets_sel = j[j_ok]
        jsel_pt  = j_pt  [j_ok]
        jsel_m   = j_mass[j_ok]

        # b-jet identification
        is_bjet  = (abs(jets_sel.eta) < 2.4) & (jets_sel.btagDeepFlavB > bWP)
        has_bjet = ak.any(is_bjet, axis=1)

        # Light-flavor jets only (no b-tag)
        lf_jets  = jets_sel[~is_bjet]
        lf_pt    = jsel_pt [~is_bjet]
        lf_mass  = jsel_m  [~is_bjet]

        # ── 9. VBF veto (leading 2 of ALL selected jets) ─────────────
        sort_all  = ak.argsort(jsel_pt, ascending=False)
        all_sorted_pt   = jsel_pt      [sort_all]
        all_sorted_eta  = jets_sel.eta [sort_all]
        all_sorted_phi  = jets_sel.phi [sort_all]
        all_sorted_mass = jsel_m       [sort_all]
        nall = ak.num(jsel_pt)

        # Pad to at least 2; fill with sentinel so cuts see the pad
        def _lead2(arr, pad):
            return ak.fill_none(ak.pad_none(arr, 2, clip=True), pad)

        j1v_pt   = ak.flatten(_lead2(all_sorted_pt,   0.0)[:, 0:1], axis=1)
        j1v_eta  = ak.flatten(_lead2(all_sorted_eta,  0.0)[:, 0:1], axis=1)
        j1v_phi  = ak.flatten(_lead2(all_sorted_phi,  0.0)[:, 0:1], axis=1)
        j1v_mass = ak.flatten(_lead2(all_sorted_mass, 0.0)[:, 0:1], axis=1)
        j2v_pt   = ak.flatten(_lead2(all_sorted_pt,   0.0)[:, 1:2], axis=1)
        j2v_eta  = ak.flatten(_lead2(all_sorted_eta,  0.0)[:, 1:2], axis=1)
        j2v_phi  = ak.flatten(_lead2(all_sorted_phi,  0.0)[:, 1:2], axis=1)
        j2v_mass = ak.flatten(_lead2(all_sorted_mass, 0.0)[:, 1:2], axis=1)

        deta_veto = abs(j1v_eta - j2v_eta)
        mjj_veto  = self._dijet_mass(j1v_pt, j1v_eta, j1v_phi, j1v_mass,
                                     j2v_pt, j2v_eta, j2v_phi, j2v_mass)
        vbf_veto_fails = (nall < 2) | (deta_veto <= 2.2) | (mjj_veto <= 450)

        # ── 10. VBF selection (light-flavor jets only) ────────────────
        sort_lf   = ak.argsort(lf_pt, ascending=False)
        lf_sorted_pt   = lf_pt        [sort_lf]
        lf_sorted_eta  = lf_jets.eta  [sort_lf]
        lf_sorted_phi  = lf_jets.phi  [sort_lf]
        lf_sorted_mass = lf_mass      [sort_lf]
        nlf = ak.num(lf_pt)

        q1_pt   = ak.flatten(_lead2(lf_sorted_pt,   0.0)[:, 0:1], axis=1)
        q1_eta  = ak.flatten(_lead2(lf_sorted_eta,  0.0)[:, 0:1], axis=1)
        q1_phi  = ak.flatten(_lead2(lf_sorted_phi,  0.0)[:, 0:1], axis=1)
        q1_mass = ak.flatten(_lead2(lf_sorted_mass, 0.0)[:, 0:1], axis=1)
        q2_pt   = ak.flatten(_lead2(lf_sorted_pt,   0.0)[:, 1:2], axis=1)
        q2_eta  = ak.flatten(_lead2(lf_sorted_eta,  0.0)[:, 1:2], axis=1)
        q2_phi  = ak.flatten(_lead2(lf_sorted_phi,  0.0)[:, 1:2], axis=1)
        q2_mass = ak.flatten(_lead2(lf_sorted_mass, 0.0)[:, 1:2], axis=1)

        deta_sel = abs(q1_eta - q2_eta)
        mjj_sel  = self._dijet_mass(q1_pt, q1_eta, q1_phi, q1_mass,
                                    q2_pt, q2_eta, q2_phi, q2_mass)
        vbf_sel_ok = (nlf >= 2) & (deta_sel > 2.2) & (mjj_sel > 450)


        pass_vbf = (~has_bjet) & (~vbf_veto_fails) & vbf_sel_ok

        # Apply VBF cuts (filter all parallel arrays)
        def _cut(arr): return arr[pass_vbf]

        events    = _cut(events)
        trig_prio = _cut(trig_prio)
        fj_pt     = _cut(fj_pt)
        fj_eta    = _cut(fj_eta)
        fj_phi    = _cut(fj_phi)
        fj_mass   = _cut(fj_mass)
        x4b       = _cut(x4b)
        pnet_v2b  = _cut(pnet_v2b)
        pnet_34   = _cut(pnet_34)
        deta_sel  = _cut(deta_sel)
        mjj_sel   = _cut(mjj_sel)
        q1_pt     = _cut(q1_pt);   q1_eta = _cut(q1_eta)
        q1_phi    = _cut(q1_phi);  q1_mass= _cut(q1_mass)
        q2_pt     = _cut(q2_pt);   q2_eta = _cut(q2_eta)
        q2_phi    = _cut(q2_phi);  q2_mass= _cut(q2_mass)

        self.cnt["pass_VBF"] += len(events)
        if len(events) == 0:
            return {}

        # ── 11. VBF sub-category ─────────────────────────────────────
        is_hi    = (deta_sel > 3.0) & (mjj_sel > 900)
        is_hi_pt = fj_pt >= 400

        vbf_cat_idx = ak.where(
            is_hi & is_hi_pt,        np.int32(3),   # VBFjjHiPtHi
            ak.where(
                is_hi & (~is_hi_pt), np.int32(2),   # VBFjjHiPtLo
                ak.where(
                    (~is_hi) & is_hi_pt, np.int32(1),  # VBFjjLoPtHi
                    np.int32(0)  # VBFjjLoPtLo
                )
            )
        )

        # ── 12. X4b regions ──────────────────────────────────────────
        is_SR      = x4b > 0.96
        is_SB      = x4b > 0.84
        is_SB_only = (x4b > 0.84) & (x4b <= 0.96)

        # ── 13. Fill output buffers ───────────────────────────────────
        def to_np(arr):
            return np.asarray(ak.to_numpy(ak.values_astype(arr, type(arr[0]) if len(arr) else float)))

        run_np  = np.asarray(ak.to_numpy(events.run))
        ls_np   = np.asarray(ak.to_numpy(events.luminosityBlock))
        evt_np  = np.asarray(ak.to_numpy(events.event))
        fpt_np  = np.asarray(ak.to_numpy(fj_pt))
        feta_np = np.asarray(ak.to_numpy(fj_eta))
        fphi_np = np.asarray(ak.to_numpy(fj_phi))
        fmas_np = np.asarray(ak.to_numpy(fj_mass))
        q1pt_np = np.asarray(ak.to_numpy(q1_pt))
        q1et_np = np.asarray(ak.to_numpy(q1_eta))
        q1ph_np = np.asarray(ak.to_numpy(q1_phi))
        q1ms_np = np.asarray(ak.to_numpy(q1_mass))
        q2pt_np = np.asarray(ak.to_numpy(q2_pt))
        q2et_np = np.asarray(ak.to_numpy(q2_eta))
        q2ph_np = np.asarray(ak.to_numpy(q2_phi))
        q2ms_np = np.asarray(ak.to_numpy(q2_mass))
        det_np  = np.asarray(ak.to_numpy(deta_sel))
        mjj_np  = np.asarray(ak.to_numpy(mjj_sel))
        x4b_np  = np.asarray(ak.to_numpy(x4b))
        pv2b_np = np.asarray(ak.to_numpy(pnet_v2b))
        p34_np  = np.asarray(ak.to_numpy(pnet_34))
        prio_np = np.asarray(ak.to_numpy(trig_prio)).astype(np.int32)
        tf_np   = np.ones(len(events), dtype=np.int32)
        cat_np  = np.asarray(ak.to_numpy(vbf_cat_idx))
        sr_np   = np.asarray(ak.to_numpy(is_SR))
        sb_np   = np.asarray(ak.to_numpy(is_SB))
        sbo_np  = np.asarray(ak.to_numpy(is_SB_only))

        for region_name, rm in [
            ("X4bSR",      sr_np.astype(bool)),
            ("X4bSB",      sb_np.astype(bool)),
            ("X4bSB_only", sbo_np.astype(bool)),
        ]:
            for cat_idx, cat_name in enumerate(VBF_CATEGORIES):
                sel = rm & (cat_np == cat_idx)
                if not np.any(sel):
                    continue
                key = f"{cat_name}_{region_name}"
                buf = self._out[key]
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
        print(f"\n{'='*70}")
        print(f"Summary for {dataset_name}:")
        print(f"{'='*70}")
        print(f"Total events processed : {self.cnt['total']}")
        print(f"Pass golden JSON       : {self.cnt['pass_golden']}")
        print(f"Pass trigger           : {self.cnt['pass_trigger']}")
        print(f"Pass fat-jet presel    : {self.cnt['pass_fatH']}")
        print(f"Pass VBF selection     : {self.cnt['pass_VBF']}")
        print("\nTree counts:")
        for cat in VBF_CATEGORIES:
            print(f"\n  {cat}:")
            for sel in X4B_SELECTIONS:
                key = f"{cat}_{sel}"
                n = len(self._out[key]["run"])
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
    print(f"Total events processed : {merged_cnt['total']:,}")
    print(f"Pass golden JSON       : {merged_cnt['pass_golden']:,}")
    print(f"Pass trigger           : {merged_cnt['pass_trigger']:,}")
    print(f"Pass fat-jet presel    : {merged_cnt['pass_fatH']:,}")
    print(f"Pass VBF selection     : {merged_cnt['pass_VBF']:,}")
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
