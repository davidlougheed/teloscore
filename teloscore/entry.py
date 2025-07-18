import argparse

from . import __version__
from .compare import Scoring1, Scoring2, compare_samples
from .plot import plot_versus


def cmd_compare(args):
    scoring = Scoring2() if args.scoring == "2" else Scoring1()
    compare_samples(scoring, args.f1, args.f2, args.out)


def cmd_plot(args):
    plot_versus(args.file)


def main():
    parser = argparse.ArgumentParser(
        prog="teloscore",
        description="A telomeric allele comparison-scoring tool for output from Telogator2.",
    )

    parser.add_argument("--version", action="version", version=__version__)

    subparsers = parser.add_subparsers()

    # Compare sub-parser -----------------------------------------------------------------------------------------------

    sp_compare = subparsers.add_parser("compare")

    default_scoring = "2"
    sp_compare.add_argument(
        "--scoring",
        choices=("1", "2"),
        help=f"Scoring system to use (1 or 2). Default: {default_scoring}",
        default=default_scoring,
    )

    sp_compare.add_argument("f1", type=str, help="First Telogator2 TSV file.")
    sp_compare.add_argument("f2", type=str, help="Second Telogator2 TSV file.")
    sp_compare.add_argument(
        "out", type=str, help="Output file to write to. This will overwrite any existing file at this path!"
    )

    sp_compare.set_defaults(func=cmd_compare)

    # Plot sub-parser --------------------------------------------------------------------------------------------------

    sp_plot = subparsers.add_parser("plot")
    sp_plot.add_argument("file", type=str, help="Output from teloscore compare function.")
    sp_plot.set_defaults(func=cmd_plot)

    # ------------------------------------------------------------------------------------------------------------------

    args = parser.parse_args()

    func = getattr(args, "func", None)

    if func is None:
        parser.parse_args(("--help",))  # will exit

    func(args)


if __name__ == "__main__":
    main()
