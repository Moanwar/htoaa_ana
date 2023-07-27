


if __name__ == "__main__":
    
    print("cal_HTRewght_QCDbGen.py:: main: {}".format(sys.argv)); sys.stdout.flush()

    parser = argparse.ArgumentParser(description='cal_HTRewght_QCDbGen')
    parser.add_argument('-era', dest='era',   type=str, default=Era_2018,                    choices=[Era_2016, Era_2017, Era_2018], required=False)
    args=parser.parse_args()
    print("args: {}".format(args))

    era              = args.era


    nRebin = 10
    sFIn_dict = {
        Era_2016: "",
        Era_2017: "",        
        Era_2018: {
            "Data": {
                "fileName": "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/UltraLegacy/PileupHistogram-goldenJSON-13tev-2018-69200ub-99bins.root",
                "histogramName": "pileup"
            },

            "MC"  : {
                "fileName": "/eos/cms/store/user/ssawant/htoaa/analysis/20230720_DataVsMC_woBtagWgt_wDeepTagbbVsL/2018/analyze_htoaa_stage1.root",
                "histogramName": "evt/$SAMPLECATEGORY/hPileup_nTrueInt_central"
            }
        }
    }

    sFOut_dict = {
        Era_2016: "",
        Era_2017: "",        
        Era_2018: "../data/correction/mc/HTSamplesStitch/HTSamplesStitchSF_2018.root"
    }

    sHistogramNameForSF = {
        "QCD_bGen": {
            "N": "evt/QCD_Incl_PSWeight/hGenLHE_HT_SelQCDbGen_central", # SF numerator histogram
            "D": "evt/QCD_bGen/hGenLHE_HT_SelQCDbGen_central", # SF denominator histogram
            "Range": [100, 1800], # axis (HT) range for SF computation
            "HTBins": [
                [ 100,  200],
                [ 200,  300],
                [ 300,  500],
                [ 500,  700],
                [ 700, 1000],
                [1000, 1500],
                [1500, 2000],
                [2000, 3000],
            ],
        }
    }

    
    sFIn  = sFIn_dict[era]
    sFOut = sFOut_dict[era]


        
    # Read input file --------------------------------
    print(f"Input file: {sFIn}")
    fIn = R.TFile(sFIn)
    if not fIn.IsOpen():
        print(f"Could not open input file: {sFIn}")
        exit(0)


    hTmp = R.TH1D("hTMp", "", 1,0,1)
    hTmp.SetDefaultSumw2()
    
    fOut = R.TFile(sFOut, 'RECREATE')
