import sys
import json
#import uproot
import uproot3 as uproot

from htoaa_Settings import * 


def calculate_lumiScale(luminosity, crossSection, sumEvents):
    lumiScale = 1
    # as crosssection is in pb and luminosity in fb
    pb_to_fb_conversionFactor = 1000
    
    if sumEvents != 0: lumiScale = luminosity * crossSection * pb_to_fb_conversionFactor / sumEvents
    return lumiScale


def setXRootDRedirector(fileName):
    if not fileName.startswith("/store/"):
        return fileName
    
    redirector_toUse = None
    for redirector in xrootd_redirectorNames:
        print(f"setXRootDRedirector():: Checking {redirector + fileName}"); sys.stdout.flush()
        '''
        with uproot.open(redirector + fileName) as file1:
            #print(f"\n{redirector + fileName}: file1.keys(): {file1.keys()}")
            #print(f"\n{redirector + fileName}: file1.keys(): {file1['Events'].numentries}"); sys.stdout.flush()

            #if file1['Events'].num_entries > 0:
            if file1['Events'].numentries > 0:
                print(f"{redirector + fileName}: file1.keys(): {file1['Events'].numentries}"); sys.stdout.flush()
                redirector_toUse = redirector
                break
        '''
        file1 = None
        try:
            file1 = uproot.open(redirector + fileName)
        except:
            print(f"setXRootDRedirector():: File open {redirector + fileName} failed"); sys.stdout.flush()
        else:
            #print(f"\n{redirector + fileName}: file1.keys(): {file1.keys()}")
            #print(f"\n{redirector + fileName}: file1.keys(): {file1['Events'].numentries}"); sys.stdout.flush()

            nEntries = file1['Events'].numentries
            file1.close()
            #if file1['Events'].num_entries > 0:
            if nEntries > 0:
                print(f"{redirector + fileName}: {nEntries}"); sys.stdout.flush()
                redirector_toUse = redirector
                break
            
    #print(f"redirector_toUse: {redirector_toUse}")
    
    return redirector_toUse + fileName


def GetDictFromJsonFile(filePath):
    # Lines starting with '#' are not read out, and also content between '/* .... */' are not read.
    # Content between " '''   ....  ''' " are not read
    # Source: https://stackoverflow.com/questions/29959191/how-to-parse-json-file-with-c-style-comments
    
    contents = ""
    fh = open(filePath)
    for line in fh:
        #cleanedLine = line.split("//", 1)[0]
        cleanedLine = line.split("#", 1)[0]
        if len(cleanedLine) > 0 and line.endswith("\n") and "\n" not in cleanedLine:
            cleanedLine += "\n"
        contents += cleanedLine
    fh.close
    
    #while "/*" in contents:
    #    preComment, postComment = contents.split("/*", 1)
    #    contents = preComment + postComment.split("*/", 1)[1]
    while "'''" in contents:
        preComment, postComment = contents.split("'''", 1)
        contents = preComment + postComment.split("'''", 1)[1]

    dictionary =  json.loads( contents )
    return dictionary


def DfColLabel_convert_bytes_to_string(df):
    cols_rename = {}
    for col in df.columns:
        if isinstance(col, (bytes, bytearray)):
            cols_rename[col] = col.decode()
    print("DfColLabel_convert_bytes_to_string:: cols_rename: {}".format(cols_rename))
    df.rename(columns=cols_rename, inplace=True)
    return df


def cut_ObjectMultiplicity(nObjects, nObjects_min=None, nObjects_max=None):
    '''
    Check if "nObjects are with nObjects_min and nObjects_max", if they are specified.
    If either of nObjects_min and nObjects_max are not specified, then corresponding condition is not checked.

    Return:
        True: if nObjects passes the condition.
        False: if nObjects fails the condition
    '''
    mask = mask_low = mask_up  = None
    if nObjects_min is not None: mask_low = (nObjects >= nObjects_min)    
    if nObjects_max is not None: mask_up  = (nObjects <= nObjects_max)

    if (nObjects_min is not None) and (nObjects_max is not None):
        mask = (mask_low and mask_up)
    elif (nObjects_min is not None):
        mask = mask_low
    else:
        mask = mask_up

    return mask



def cut_ObjectPt(objects_Pt, PtThrsh_Lead=None, PtThrsh_Sublead=None, PtThrsh_Third=None, PtThrsh_Fourth=None, PtThrsh_Fifth=None):
    '''
    Check Objects Pt is above their resepctive thresholds.
    If PtThrsh_<rank> is not set, then their Pt condition is not checked.

    Return:
        True: All objects' Pt is about respective threshold
        False: Else false
    '''
    print("objects_Pt ({}) : {}".format(type(objects_Pt),  objects_Pt))
    condition = True
    if                                    objects_Pt[0] < PtThrsh_Lead:     condition = False    
    if (PtThrsh_Sublead is not None) and (objects_Pt[1] < PtThrsh_Sublead): condition = False
    if (PtThrsh_Third   is not None) and (objects_Pt[2] < PtThrsh_Third):   condition = False
    if (PtThrsh_Fourth  is not None) and (objects_Pt[3] < PtThrsh_Fourth):  condition = False
    if (PtThrsh_Fifth   is not None) and (objects_Pt[4] < PtThrsh_Fifth):   condition = False

    return condition
    
def cut_ObjectPt_1(objects_Pt, PtThrshs):
    print("objects_Pt ({}) : {}, \t\t PtThrshs ({}) : {}".format(type(objects_Pt),  objects_Pt, type(PtThrshs), PtThrshs))

    return True

def cut_ObjectEta(objects_Eta, EtaThrsh, nObjects):
    '''
    Check Objects abs(Eta) is greater than thresholds set in EtaThrshs list.

    Return:
        True: All objects' Eta is about respective threshold
        False: Else false
    '''    
    condition = True
    for iObject in range(nObjects):
        if abs(objects_Eta[iObject]) > EtaThrsh:
            condition = False
            break
    return condition
