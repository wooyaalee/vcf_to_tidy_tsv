import sys
from collections import defaultdict
import gzip
import argparse

base_header = "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tTUMOR\tNORMAL\tCALLER"

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--tumor", type=str, dest="tumor", help = "Tumor Sample Label", default=None)
    parser.add_argument("-n", "--normal", type=str, dest="normal", help = "Normal Sample Label", default=None)
    parser.add_argument("-c", "--caller", type = str, dest="caller", help = "The caller to label calls with.", required=True)
    parser.add_argument("-v", "--vcf", type = str, dest="vcf", help = "A variant call file to transform.", required=True)
    
    ## parser.add_argument("-a", "--annotations",
    ## type=str, dest="annotations",
    ## help="A file containing annotations to apply to ALL VARIANTS, where each line in the file is KEY:VALUE")

    return parser.parse_args()

## #CHROM POS REF ALT TUMOR NORMAL CALLER FORMAT:<key> ... INFO:<key> ...
def print_header(header_d): 
    h = base_header + "\t" + "\t".join(sorted([header_d[x] for x in header_d]))
    print( "#" + h)
    return h

def make_header_index_d(header):
    d = defaultdict(int)
    inverted_d = defaultdict(str)
    index = 0
    for h in header.strip().strip("#").split("\t"):
        d[h] = index
        inverted_d[index] = h
        index += 1
    return d, inverted_d


if __name__ == "__main__":
    
## We need to tweak the way our header works a bit
## We're going to set the column names to 
## INFO:<key>, FORMAT:<key>, etc
    header_d = defaultdict(str)
    header_index_d = defaultdict(int)
    info_flags = defaultdict(str)
    
    args = parse_args()
    
    callerLabel = args.caller
    pair = args.vcf.strip("./").split(".")[0]
    tumorLabel = pair
    if args.tumor is not None:
        tumorLabel = args.tumor
    normalLabel = pair
    if args.normal is not None:
        normalLabel = args.normal


    headerTripped = False

    ifi = None

    if ".gz" in args.vcf:
        ifi = gzip.open(args.vcf, "r")
    else:
        ifi = open(args.vcf, "r")
    for line in ifi:
        if line.startswith("##"):
            if "##INFO" in line:
                idVal = [i for i in line.replace("##INFO=", "").strip("<>").split(",") if "ID" in i][0].split("=")[1]
                if "Type=Flag" in line:
                    info_flags["INFO" + idVal] = ""
                header_d["INFO:" + idVal] = "INFO:" + idVal
            elif "##FORMAT" in line:
                idVal = [i for i in line.replace("##FORMAT=", "").strip("<>").split(",") if "ID" in i][0].split("=")[1]
                header_d["FORMAT:" + idVal] = "FORMAT:" + idVal
        else:
            if headerTripped:
                tokens = line.strip().split("\t")
                ## This dictionary holds our line outputs
                ##  in sorted order (because the index is the output index)
                ## Remember, our base header is:
                ## "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tTUMOR\tNORMAL\tCALLER\t"
                ## Get our first few fields, which are fixed, and also add our caller / tumor / normal info
                outputs = defaultdict(str)
                outputs[0] = tokens[0]
                outputs[1] = tokens[1]
                outputs[2] = tokens[2]
                outputs[3] = tokens[3]
                outputs[4] = tokens[4]
                outputs[5] = tokens[5]
                outputs[6] = tokens[6]
                outputs[7] = tumorLabel
                outputs[8] = normalLabel
                outputs[9] = args.caller
                ## break out our info fields
                ## These have KEY=VALUE style, so easier to grab them
                infos = tokens[7].strip().split(";")
                ## Remember 1 special case:
                ## handle flag vals, which, when split by "=",
                ## are only one element, and may not be present.
                ## We label these like booleans.
                ## everything else is a stringifiable value
                flags = defaultdict(bool)
                for i in infos:
                    splits = i.split("=")
                    rawKey = splits[0]
                    key = ":".join(["INFO", rawKey])
                    if len(splits) > 1:
                        val = splits[1]
                        outputs[header_index_d[key]] = val
                    else:
                        flags[key] = True
                for i in info_flags:
                    if key in flags:
                        outputs[header_index_d[key]] = "TRUE"
                    else:
                        outputs[header_index_d[key]] = "FALSE"

                ## Fill missing values
                for i in inverted_header_index_d:
                    if i not in outputs:
                        outputs[i] = "NA"


                ## Lastly, print our new TSV style line
                # print("HEADER", header_d)
                # print("INDEX", header_index_d)
                # print("INVERTED", inverted_header_index_d)
                # print("OUTPUTS", outputs)
                # print("TOKENS", tokens)
                # exit(1)
                print("\t".join(outputs[i] for i in sorted(outputs)))
            else:
                header = print_header(header_d)
                header_index_d, inverted_header_index_d = make_header_index_d(header)
                headerTripped = True

    ifi.close()


