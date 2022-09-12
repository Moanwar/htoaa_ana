import json




def calculate_lumiScale(luminosity, crossSection, sumEvents):
    lumiScale = 1
    if sumEvents != 0: lumiScale = luminosity * crossSection / sumEvents
    return lumiScale


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
