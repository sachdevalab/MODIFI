#!/usr/bin/env python
import argparse

import pandas as pd

# import polars as pl
# import skbio.io
import smart_open
from Bio.SeqIO.FastaIO import SimpleFastaParser

# import csv


# def filter_blast(blastn_in, spacer_in, max_mismatch, blast_out):
#     spacer_id_to_seq = {}

#     for header, seq in SimpleFastaParser(smart_open.open(spacer_in)):
#         id_ = header.split()[0]

#         spacer_id_to_seq[id_] = seq

#     header = (
#         "qseqid",
#         "sseqid",
#         "pident",
#         "length",
#         "mismatch",
#         "gapopen",
#         "qstart",
#         "qend",
#         "sstart",
#         "send",
#         "evalue",
#         "bitscore",
#         "qlen",
#         "nident",
#     )

#     df_blast = pl.read_csv(
#         blastn_in, separator="\t", has_header=False, new_columns=header
#     )

#     df_blast = df_blast.with_columns(
#         pl.col("qseqid")
#         .map_elements(lambda x: "_".join(x.split("_CRISPR_")[:-1]))
#         .alias("query_contig_id")
#     )
#     df_blast = df_blast.filter(pl.col("query_contig_id") != pl.col("sseqid"))
#     df_blast = df_blast.with_columns(
#         (pl.col("qlen") - pl.col("nident")).alias("full_mismatch")
#     )
#     df_blast = df_blast.filter(pl.col("full_mismatch") <= max_mismatch)
#     df_blast = df_blast.with_columns(
#         (pl.col("nident") / pl.col("qlen") * 100).alias("full_percent_id")
#     )
#     df_blast = df_blast.with_columns(
#         pl.col("qseqid").map_dict(spacer_id_to_seq).alias("spacer_seq")
#     )

#     # df_blast = df_blast.collect(streaming=True)

#     df_blast = df_blast.select(
#         [
#             "qseqid",
#             "query_contig_id",
#             "sseqid",
#             "full_mismatch",
#             "full_percent_id",
#             "spacer_seq",
#         ]
#     )

#     df_blast = df_blast.rename(
#         {
#             "qseqid": "query_spacer_id",
#             "sseqid": "target_id",
#         }
#     )

#     # print("spacer.test.tsv")

#     return df_blast


def filter_blast(blast_in, spacer_in, max_mismatch, blast_out):
    spacer_id_to_seq = {}

    for header, seq in SimpleFastaParser(smart_open.open(spacer_in)):
        id_ = header.split()[0]

        spacer_id_to_seq[id_] = seq

    # blast_df = skbio.io.read(smart_open.open(blast_in), into=pd.DataFrame)
    #

    header = (
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "qlen",
        "nident",
    )

    blast_df = pd.read_csv(blast_in, sep="\t", header=None, names=header)

    blast_df["qseqid"] = blast_df["qseqid"].astype(str)
    blast_df["sseqid"] = blast_df["sseqid"].astype(str)

    blast_df["query_contig_id"] = blast_df["qseqid"].apply(
        lambda x: "_".join(x.split("_CRISPR_")[:-1])
    )

    blast_df = blast_df[blast_df["query_contig_id"] != blast_df["sseqid"]]

    blast_df["full_mismatch"] = blast_df["qlen"] - blast_df["nident"]
    blast_df = blast_df[blast_df["full_mismatch"] <= max_mismatch]
    blast_df["full_percent_id"] = blast_df["nident"] / blast_df["qlen"] * 100

    blast_df["spacer_seq"] = blast_df["qseqid"].apply(lambda x: spacer_id_to_seq[x])

    blast_df = blast_df[
        [
            "qseqid",
            "query_contig_id",
            "sseqid",
            "full_mismatch",
            "full_percent_id",
            "spacer_seq",
        ]
    ]

    blast_df = blast_df.rename(
        columns={
            "qseqid": "query_spacer_id",
            "sseqid": "target_id",
        }
    )

    return blast_df


# def filter_blast(blast_in, spacer_in, max_mismatch, blast_out):
#     spacer_id_to_seq = {}

#     for record in SeqIO.parse(smart_open.open(spacer_in), "fasta"):
#         spacer_id_to_seq[record.id] = str(record.seq)

#     blast_df = skbio.io.read(smart_open.open(blast_in), into=pd.DataFrame)
#     #
#     blast_df["query_contig_id"] = blast_df["qaccver"].apply(
#         lambda x: "_".join(x.split("_CRISPR_")[:-1])
#     )

#     blast_df = blast_df[blast_df["query_contig_id"] != blast_df["saccver"]]

#     blast_df["full_mismatch"] = blast_df["qlen"] - blast_df["nident"]
#     blast_df = blast_df[blast_df["full_mismatch"] <= max_mismatch]
#     blast_df["full_percent_id"] = blast_df["nident"] / blast_df["qlen"] * 100

#     blast_df["spacer_seq"] = blast_df["qaccver"].apply(lambda x: spacer_id_to_seq[x])

#     blast_df = blast_df[
#         [
#             "qaccver",
#             "query_contig_id",
#             "saccver",
#             "full_mismatch",
#             "full_percent_id",
#             "spacer_seq",
#         ]
#     ]

#     blast_df = blast_df.rename(
#         columns={
#             "qaccver": "query_spacer_id",
#             "saccver": "target_id",
#         }
#     )

#     return blast_df


# def main():
#     df_blast = filter_blast(
#         args.blast_in, args.spacer_in, args.max_mismatch, args.blast_out
#     )

#     df_blast = df_blast.write_csv(args.blast_out, separator="\t")


def main():
    blast_df = filter_blast(
        args.blast_in, args.spacer_in, args.max_mismatch, args.blast_out
    )

    blast_df = blast_df.to_csv(args.blast_out, sep="\t", index=False)


if __name__ == "__main__":
    argparser = argparse.ArgumentParser()

    argparser.add_argument(
        "-i",
        "--in",
        action="store",
        dest="blast_in",
        required=True,
    ),

    argparser.add_argument(
        "-o",
        "--out",
        action="store",
        dest="blast_out",
        required=True,
    ),

    argparser.add_argument(
        "-s",
        "--spacers",
        action="store",
        dest="spacer_in",
        required=True,
    ),

    argparser.add_argument(
        "-m",
        "--max_mismatch",
        action="store",
        dest="max_mismatch",
        required=False,
        default=1,
        type=int,
    )

    args = argparser.parse_args()

    main()
