#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import Bio
from Bio import SeqIO
import argparse


def parse_args():
    """ 
	Gets the arguments from the command line.
	"""

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--fasta", required=True, help="""Fasta file to filter."""
    )
    parser.add_argument(
        "-o",
        "--output",
        required=False,
        default="output.fasta",
        help="""Output fasta [default: output.fasta]""",
    )
    parser.add_argument("-m", "--min", required=True, help="""Minimum read length.""")

    return parser.parse_args()


def main():
    args = parse_args()
    fasta_file = args.fasta
    output_file = args.output
    min_len = args.min

    record = SeqIO.parse(fasta_file, "fasta")
    output = open(output_file, "w")

    for sequence in record:
        if len(sequence.seq) >= min_len:
            output.write(">{0}\n{1}\n".format(sequence.id, sequence.seq))


if __name__ == "__main__":
    main()
