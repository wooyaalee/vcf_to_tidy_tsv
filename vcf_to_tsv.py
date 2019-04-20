import sys
from collections import defaultdict
import gzip
import argparse

base_header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tTUMOR\tNORMAL\tCALLER\t"

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--tumor", type=str, dest="tumor", help = "Tumor Sample Label", default=None)
    parser.add_argument("-n", "--normal", type=str, dest="normal", help = "Normal Sample Label", default=None)
    parser.add_argument("-c", "--caller", type = str, dest="caller", help = "The caller to label calls with.", required=True)
    parser.add_argument("-v", "--vcf", type = str, dest="vcf", help = "A variant call file to transform.", required=True)

    return parser.parse_args()

## #CHROM POS REF ALT TUMOR NORMAL CALLER FORMAT:<key> ... INFO:<key> ...
def print_header(header_d): 
    h = base_header + "\t".join(sorted([header_d[x] for x in header_d]))
    print(h)
    return h

def make_header_index_d(header):
    d = defaultdict(int)
    index = 0
    for h in header.strip().split("\t"):
        d[h] = index
        index += 1
    return d


if __name__ == "__main__":
    
## We need to tweak the way our header works a bit
## We're going to set the column names to 
## INFO:<key>, FORMAT:<key>, etc
    header_d = defaultdict(str)
    header_index_d = defaultdict(int)
    
    args = parse_args()
    
    callerLabel = args.caller
    pair = args.vcf.split(".")[0]
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
                header_d[idVal] = "INFO:" + idVal
            elif "##FORMAT" in line:
                idVal = [i for i in line.replace("##FORMAT=", "").strip("<>").split(",") if "ID" in i][0].split("=")[1]
                header_d[idVal] = "FORMAT:" + idVal
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
                ## Special case: handle flag vals, which, when split by "=",
                ## are only one element



                ## Lastly, print our new TSV style line
                print("\t".join(outputs[i] for i in sorted(outputs)))
            else:
                header = print_header(header_d)
                header_index_d = make_header_index_d(header)
                headerTripped = True

    ifi.close()


