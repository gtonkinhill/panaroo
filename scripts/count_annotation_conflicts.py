import argparse
import sys, os
import re
# import Levenshtein as lev


def read_gene_data(inputfile):

    gene_data = {}
    with open(inputfile, 'r') as infile:
        next(infile)
        for line in infile:
            line = line.strip().split(",")
            gene_data[line[3]] = (len(line[4]), line[-2])
    return (gene_data)


def count_differences(gene_data,
                      pa_file,
                      lendiff=0.8,
                      col_skip=0,
                      sep=",",
                      method="panaroo"):
    anno_conflicts = 0
    record_anno_conflicts = []
    with open(pa_file, 'r') as infile:
        next(infile)
        for line in infile:
            realline = line
            if method == "roary":
                line = re.split('","', line.strip())[col_skip:]
            else:
                line = line.strip().split(sep)[col_skip:]

            cluster = []
            for l in line:
                if method == "pirate":
                    if ":" in l: continue
                elif method == "panaroo":
                    if ";" in l: continue

                l = re.split("[;:\t]", l)
                for g in l:
                    g = g.strip("(").strip(")").strip('"')
                    if "refound" in g: continue
                    if g == "": continue
                    if g not in gene_data:
                        print(line)
                        print(realline)
                    cluster.append(g)

            max_len = -1
            for g in cluster:
                max_len = max(max_len, gene_data[g][0])

            annotations = set()
            for g in cluster:
                # if gene_data[g][0] <= (lendiff*max_len):
                #     continue
                if gene_data[g][1] in ["", "hypothetical protein"]:
                    continue
                annotations.add(gene_data[g][1])

            annotations = set(
                [re.sub('[0-9|_]', '', i.lower()) for i in annotations])
            if len(annotations) > 1:
                record_anno_conflicts.append(list(annotations))
                anno_conflicts += len(annotations) - 1

            # annotations = sorted(list(annotations))
            # if len(annotations)>1:
            #     for a in annotations[1:]:
            #         if lev.ratio(annotations[0].lower(), a.lower()) < lev_thresh:
            #             anno_conflicts += 1
            #             record_anno_conflicts.append([annotations[0].lower(), a.lower()])

    for a in record_anno_conflicts:
        print(a)

    return anno_conflicts


def main():
    parser = argparse.ArgumentParser(
        description="""Counts annotation conflicts.""")

    parser.add_argument("-g",
                        "--gene_data",
                        dest="gene_data",
                        required=True,
                        help="gene data file output by Panaroo",
                        type=str)

    parser.add_argument("-p"
                        "--pa",
                        dest="pa_file",
                        help="Presence absence file",
                        required=True)

    parser.add_argument(
        "--method",
        dest="method",
        help="Algorithm used to produce p/a file for formatting",
        type=str,
        choices=['panaroo', 'roary', 'pirate', 'panx', 'cogsoft'],
        default="panaroo")

    args = parser.parse_args()

    gene_data = read_gene_data(args.gene_data)

    if args.method == "pirate":
        col_skip = 22
        sep = "\t"
    elif args.method == "panaroo":
        col_skip = 14
        sep = ","
    elif args.method == "roary":
        col_skip = 14
        sep = ","
    elif args.method == "panx":
        col_skip = 0
        sep = ","
    elif args.method == "cogsoft":
        col_skip = 0
        sep = ","

    anno_conlfict_count = count_differences(gene_data,
                                            args.pa_file,
                                            col_skip=col_skip,
                                            sep=sep,
                                            method=args.method)

    print("Conflicting annotations: ", anno_conlfict_count)

    return


if __name__ == '__main__':
    main()