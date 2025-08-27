#!/usr/bin/env python3
import argparse

from Bio.SeqIO.FastaIO import SimpleFastaParser


def rename_contigs(in_assembly_path, sample_name, out_assembly_path):
    with open(out_assembly_path, "w") as out_assembly_path_fh:
        for header, seq in SimpleFastaParser(open(in_assembly_path)):
            id_ = header.split()[0]

            if (
                id_.startswith("p")
                and "ptg" in id_
                and (id_.endswith("l") or id_.endswith("c"))
            ):
                contig_num = id_.split("ptg")[1].replace("l", "").replace("c", "")
                contig_topology = id_[-1].upper()

                contig_num = str(int(contig_num))

                edit_id = "_".join((sample_name, contig_num, contig_topology))

            else:
                contig_num = id_.split("_")[1]

                edit_id = "_".join((sample_name, contig_num))

            out_fa = ">" + edit_id + "\n" + seq
            out_assembly_path_fh.write(out_fa + "\n")


def main():
    rename_contigs(args.in_assembly, args.sample, args.out_assembly)


if __name__ == "__main__":
    argparser = argparse.ArgumentParser()

    argparser.add_argument(
        "-i",
        "--in",
        action="store",
        dest="in_assembly",
        required=True,
        help="Input sequences in FASTA format",
    )

    argparser.add_argument(
        "-s",
        "--sample",
        action="store",
        dest="sample",
        required=True,
        help="Sample name",
    )

    argparser.add_argument(
        "-o",
        "--out",
        action="store",
        dest="out_assembly",
        required=True,
        help="Output sequences in FASTA format",
    )

    args = argparser.parse_args()

    main()
