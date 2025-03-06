import argparse
import os
import sys
import warnings
import logging
from Bio import BiopythonWarning
from .genome_annotator import annotate_genome


warnings.filterwarnings("ignore", category=BiopythonWarning)

logging.basicConfig(
    stream=sys.stdout,
    format="%(asctime)s - %(message)s",
    datefmt="%d-%b-%y %H:%M:%S",
    level=logging.INFO,
)


def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error(f"Input file ({arg}) not found!")

    else:
        return arg


def main():
    parser = argparse.ArgumentParser()

    #############
    # Annotator #
    #############

    parser.add_argument(
        "-p",
        "--pfam",
        metavar="\b",
        type=lambda x: is_valid_file(parser, x),
        help="HMM Profile(s)",
        required=True,
    )

    parser.add_argument(
        "-s",
        "--sequences",
        metavar="\b",
        type=lambda x: is_valid_file(parser, x),
        help="Path to genome",
        required=True,
    )

    parser.add_argument(
        "-w",
        "--window",
        metavar="\b",
        default=100000,
        type=int,
        help="Window size (Default: 100000 bp)",
        required=False,
    )

    parser.add_argument(
        "-l",
        "--overlap",
        metavar="\b",
        default=100,
        type=int,
        help="overlap between windows (Default: 100 bp)",
        required=False,
    )

    parser.add_argument(
        "-c",
        "--cores",
        metavar="\b",
        default=1,
        type=int,
        help="Number of Cores. Note! This should a 1/4th of the total cores you wish to allocate to the tool",
        required=False,
    )

    parser.add_argument(
        "-o",
        "--output",
        metavar="\b",
        type=str,
        help="Output txt file of HMM annotations",
        required=True,
    )

    parser.add_argument(
        "-b",
        "--bed",
        metavar="\b",
        default=None,
        type=str,
        help="Output BED file of HMM annotations",
        required=False,
    )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    annotate_genome(
        pfam=args.pfam,
        genome=args.sequences,
        window_size=args.window,
        overlap=args.overlap,
        cores=args.cores,
        output=args.output,
        bed=args.bed,
    )


if __name__ == "__main__":
    main()
    sys.exit(0)
