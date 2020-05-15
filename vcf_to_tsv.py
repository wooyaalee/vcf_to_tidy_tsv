import sys
from collections import defaultdict
import gzip
import argparse

base_header = "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tTUMOR\tNORMAL\tCALLER\tSAMPLE"
#base_header = "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tCALLER"

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--tumor", type=str,
            dest="tumor", help = "Tumor Sample Label", default=None)
    parser.add_argument("-n", "--normal", type=str,
            dest="normal", help = "Normal Sample Label", default=None)
    parser.add_argument("-s", "--sample", type=str,
            dest="sample", help = "A string to use in the sample field.", default=None)
    parser.add_argument("-c", "--caller", type = str,
            dest="caller", help = "The caller to label calls with.", required=False, default="CALLER")
    parser.add_argument("-v", "--vcf", type = str,
            dest="vcf", help = "A variant call file to transform.", required=True)
    
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


## Makes a dictionary of index->sample name,
## so that we can label our format fields per sample
def count_formats(chromLine):
    sample_d = defaultdict(str)
    tokens = chromLine.strip().split("\t")
    numSamples = len(tokens) - 8
    count = 0
    if numSamples == 0:
        return sample_d
    for i in tokens[9:]:
        sample_d[count] = i
        count += 1
    return sample_d


## Duplicates format columns, so that 
## they're represented once for each sample
## Their final form will look like:
## FORMAT:<key>:<sample>
def correct_header_format(header_d, sample_d):
    reheader_d = defaultdict(str)
    for h in header_d:
        if "FORMAT" in h:
            for i in sample_d:
                mod = h + ":" + sample_d[i]
                reheader_d[mod] = mod
        else:
            reheader_d[h] = h
    return reheader_d


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

    sample_d = None

    if ".gz" in args.vcf:
        ifi = gzip.open(args.vcf, "rt")
    else:
        ifi = open(args.vcf, "r")
    for line in ifi:
        if line.startswith("##"):
            if "##INFO" in line:
                idVal = [i for i in line.replace("##INFO=", "").strip("<>").split(",") if "ID" in i][0].split("=")[1]
                if "Type=Flag" in line:
                    info_flags["INFO:" + args.caller + ":" + idVal] = ""
                header_d["INFO:" + idVal] = "INFO:" + args.caller + ":" + idVal
            elif "##FORMAT" in line:
                idVal = [i for i in line.replace("##FORMAT=", "").strip("<>").split(",") if "ID" in i][0].split("=")[1]
                header_d["FORMAT:" + idVal] = "FORMAT:" + args.caller + ":" + idVal
        elif line.startswith("#CHROM"):
            sample_d = count_formats(line)
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
                    key = ":".join(["INFO", args.caller, rawKey])
                    if len(splits) > 1:
                        val = splits[1]
                        outputs[header_index_d[key]] = val
                    else:
                        flags[key] = True
                for i in info_flags:
                    if i in flags:
                        outputs[header_index_d[key]] = "TRUE"
                    else:
                        outputs[header_index_d[key]] = "FALSE"



                ## Format field time.
                ## This counter tracks the number of samples
                ## that are in the VCF, so that we can have a fighting
                ## shot at mapping them back later.
                SAMPLE_FIELD_COUNTER = 0
                format_tags = tokens[8].strip().split(":")
                if len(tokens) > 9 and sample_d is not None:
                    for i in tokens[9:]:
                        fields = i.split(":")
                        #print(format_tags)
                        #print(fields)
                        for f in range(0, len(fields)):
                            key = "FORMAT:" + args.caller + ":" + format_tags[f] + ":" + sample_d[SAMPLE_FIELD_COUNTER]
                            outputs[header_index_d[key]] = fields[f]
                        SAMPLE_FIELD_COUNTER += 1
                        
                ## Fill missing values
                for i in inverted_header_index_d:
                    if i not in outputs:
                        outputs[i] = "NA"


                # print("HEADER", header_d)
                # print("INDEX", header_index_d)
                # print("INVERTED", inverted_header_index_d)
                # print("OUTPUTS", outputs)
                # print("TOKENS", tokens)
                # exit(1)
                
                ## Lastly, print our new TSV style line
                print("\t".join(outputs[i] for i in sorted(outputs)))
            else:
                header_d = correct_header_format(header_d, sample_d)
                header = print_header(header_d)
                header_index_d, inverted_header_index_d = make_header_index_d(header)
                headerTripped = True

    ifi.close()


