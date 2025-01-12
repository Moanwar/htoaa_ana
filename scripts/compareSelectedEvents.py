import os


if __name__ == '__main__':
    sFSelEvents_0 = "yields_HtoAA_gg0lHi_AWB.txt"

    selRunLsEvent_list0 = []
    with open(sFSelEvents_0) as fSelEvents:
        iLine = 0
        for line_ in fSelEvents:
            run_    = 1
            ls_     = int(line_.split(',')[0])
            event_  = int(line_.split(',')[1])
            if iLine < 5:
                print(f"{line_ = } \t\t >> {run_ = }, {ls_ = }, {event_ = }, ")
            iLine += 1
            selRunLsEvent_list0.append( (run_, ls_, event_) )

    print(f"selRunLsEvent_list0 ({len(selRunLsEvent_list0)}): {selRunLsEvent_list0} \n ({len(selRunLsEvent_list0)})")



    sFSelEvents_1 = "selectedRunLsEvents_SUSY_GluGluH_01J_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8.txt"

    selRunLsEvent_list1 = []
    with open(sFSelEvents_1) as fSelEvents:
        iLine = 0
        for line_ in fSelEvents:
            run_    = 1
            ls_     = int(line_.split('*')[0])
            event_  = int(line_.split('*')[1])
            if iLine < 5:
                print(f"{line_ = } \t\t >> {run_ = }, {ls_ = }, {event_ = }, ")
            iLine += 1
            selRunLsEvent_list1.append( (run_, ls_, event_) )

    print(f"selRunLsEvent_list1 ({len(selRunLsEvent_list1)}): {selRunLsEvent_list1} \n ({len(selRunLsEvent_list1)})")




    selRunLsEvent_list0 = set(selRunLsEvent_list0)
    selRunLsEvent_list1 = set(selRunLsEvent_list1)
    
    print(f"\n\n\n\nComparing selected events in file0 ({sFSelEvents_0}) and file1 ({sFSelEvents_1})")
    print(f"\n\nEvents in file0 and not in file1 ({len(selRunLsEvent_list0 - selRunLsEvent_list1)}): \n{selRunLsEvent_list0 - selRunLsEvent_list1}")
    print(f"\n\nEvents not in file0 and in file1 ({len(selRunLsEvent_list1 - selRunLsEvent_list0)}): \n{selRunLsEvent_list1 - selRunLsEvent_list0}")
    
