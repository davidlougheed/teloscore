import argparse

from . import __version__
from .compare import compare_samples


def main():
    parser = argparse.ArgumentParser(
        prog="fasta-checksum-utils",
        description="A library and command-line utility for checksumming FASTA files and individual contigs.",
    )

    parser.add_argument("--version", action="version", version=__version__)

    parser.add_argument("f1", type=str, help="First Telogator2 TSV file.")
    parser.add_argument("f2", type=str, help="Second Telogator2 TSV file.")
    parser.add_argument(
        "out", type=str, help="Output file to write to. This will overwrite any existing file at this path!"
    )

    args = parser.parse_args()

    compare_samples(args.f1, args.f2, args.out)


if __name__ == "__main__":
    main()
