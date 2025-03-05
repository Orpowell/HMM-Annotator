import argparse
import os
import sys

from genome_annotator import annotate_genome

def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error(f"Input file ({arg}) not found!")

    else:
        return arg


def main():
    parser = argparse.ArgumentParser()
    sub_parsers = parser.add_subparsers(dest="command", prog="Tool")

    #############
    # Annotator #
    #############

    annotator_parser = sub_parsers.add_parser(
        "annotator", help="Annotate Pfam domains across the genome"
    )

    annotator_parser.add_argument(
        "-p",
        "--pfam",
        metavar="\b",
        type=lambda x: is_valid_file(parser, x),
        help="Pfam hmm",
        required=True,
    )

    annotator_parser.add_argument(
        "-g",
        "--genome",
        metavar="\b",
        type=lambda x: is_valid_file(parser, x),
        help="Path to genome",
        required=True,
    )

    annotator_parser.add_argument(
        "-o",
        "--output",
        metavar="\b",
        type=str,
        help="Path to output directory",
        required=True,
    )

    annotator_parser.add_argument(
        "-w",
        "--window",
        metavar="\b",
        default=100000,
        type=int,
        help="Path to output directory",
        required=False,
    )

    annotator_parser.add_argument(
        "-l",
        "--overlap",
        metavar="\b",
        default=100,
        type=int,
        help="Path to output directory",
        required=False,
    )

    annotator_parser.add_argument(
        "-c",
        "--cores",
        metavar="\b",
        default=1,
        type=int,
        help="Number of Cores. Note! This should a 1/4th of the cores you wish to allocate to the process",
        required=False,
    )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    if args.command == "annotator":
        annotate_genome(
            pfam=args.pfam,
            genome=args.genome,
            outdir=args.output,
            window_size=args.window,
            overlap=args.overlap,
            cores=args.cores,
        )


if __name__ == "__main__":
    main()
