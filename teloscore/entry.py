import cappa
import sys

from dataclasses import dataclass
from pathlib import Path
from typing import Annotated

from . import __version__
from .compare import compare_samples
from .scoring import Scoring1, Scoring2
from .plot import plot_versus


default_scoring = "2"


@dataclass
class Compare:
    f1: Annotated[Path, cappa.Arg(help="First Telogator2 TSV file.")]
    f2: Annotated[Path, cappa.Arg(help="Second Telogator2 TSV file.")]
    out: Annotated[Path, cappa.Arg(help="Output file to write to. This will overwrite any existing file at this path!")]

    scoring: Annotated[
        str,
        cappa.Arg(
            choices=["1", "2"],
            long=True,
            help=f"Scoring system to use (1 or 2). Default: {default_scoring}"
        ),
    ] = default_scoring

    def __call__(self):
        scoring = Scoring2() if self.scoring == "2" else Scoring1()
        compare_samples(scoring, self.f1, self.f2, self.out)


@dataclass
class Plot:
    file: Annotated[Path, cappa.Arg(help="Output from teloscore compare function.")]

    def __call__(self):
        plot_versus(self.file)


@dataclass
class Version:
    def __call__(self):
        print(__version__)


@cappa.command(name="teloscore", description="A telomeric allele comparison-scoring tool for output from Telogator2.")
@dataclass
class TeloScore:
    cmd: cappa.Subcommands[Compare | Plot | Version]


def cmd_compare(args):
    scoring = Scoring2() if args.scoring == "2" else Scoring1()
    compare_samples(scoring, args.f1, args.f2, args.out)


def cmd_plot(args):
    plot_versus(args.file)


def main(argv: list[str] | None = None):
    cappa.invoke(TeloScore, argv=argv or sys.argv[1:])


if __name__ == "__main__":
    main()
