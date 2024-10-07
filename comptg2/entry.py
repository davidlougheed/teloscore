import argparse

from . import __version__
from .compare import compare_samples
from .plot import plot_versus


def cmd_compare(args):
    compare_samples(args.f1, args.f2, args.out)


def cmd_plot(args):
    plot_versus(args.file)


def main():
    parser = argparse.ArgumentParser(
        prog="comptg2",
        description="A telomeric allele comparison-scoring tool for output from Telogator2.",
    )

    parser.add_argument("--version", action="version", version=__version__)

    subparsers = parser.add_subparsers()

    # Compare sub-parser -----------------------------------------------------------------------------------------------

    sp_compare = subparsers.add_parser("compare")

    sp_compare.add_argument("f1", type=str, help="First Telogator2 TSV file.")
    sp_compare.add_argument("f2", type=str, help="Second Telogator2 TSV file.")
    sp_compare.add_argument(
        "out", type=str, help="Output file to write to. This will overwrite any existing file at this path!"
    )

    sp_compare.set_defaults(func=cmd_compare)

    # Plot sub-parser --------------------------------------------------------------------------------------------------

    sp_plot = subparsers.add_parser("plot")
    sp_plot.add_argument("file", type=str, help="Output from comptg2 compare function.")
    sp_plot.set_defaults(func=cmd_plot)

    # ------------------------------------------------------------------------------------------------------------------

    args = parser.parse_args()

    func = getattr(args, "func", None)

    if func is None:
        parser.parse_args(("--help",))  # will exit

    func(args)


if __name__ == "__main__":
    main()
