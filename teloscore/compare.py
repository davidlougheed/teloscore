import csv
import logging
from io import StringIO
from pathlib import Path
from typing import Iterable, TextIO

from .scoring import BaseScoringSystem
from .types import TeloType

__all__ = [
    "BaseScoringSystem",
    "compare_samples",
]

logging.basicConfig(level=logging.INFO)

logger = logging.getLogger("teloscore")


def build_telo_from_row(scoring: BaseScoringSystem, row: dict) -> TeloType:
    tvr = row["tvr_consensus"]
    return {
        "arm": row["#chr"],
        "allele_id": row["allele_id"],
        "tvr_consensus": tvr,
        "tvr_consensus_encoded": scoring.encode_seq(tvr),
    }


def _fmt_allele(telo: TeloType) -> str:
    return f"({telo['allele_id']}) {telo['arm']}"


def _fmt_alignment(seq1: str, seq2: str, cigar: Iterable[tuple[int, int]]) -> str:
    chars1 = []
    line_chars = []
    chars2 = []

    qs, dbs = (seq1, seq2) if len(seq1) <= len(seq2) else (seq2, seq1)

    s1 = list(qs)
    s2 = list(dbs)

    for op, count in cigar:
        # CIGAR operations are detailed here: https://samtools.github.io/hts-specs/SAMv1.pdf section 1.4 list item 6
        if op in (0, 7, 8):
            for _ in range(count):
                chars1.append(s1.pop(0))
                line_chars.append("|" if op != 8 else "X")
                chars2.append(s2.pop(0))
        elif op == 1:
            for _ in range(count):
                chars1.append(s1.pop(0))
                line_chars.append(" ")
                chars2.append("-")
        elif op in (2, 3):
            for _ in range(count):
                chars1.append("-")
                line_chars.append(" ")
                chars2.append(s2.pop(0))
        elif op == 4:
            for _ in range(count):
                s1.pop(0)

    return f"{''.join(chars1)}\n{''.join(line_chars)}\n{''.join(chars2)}"


def _read_sample_files(scoring: BaseScoringSystem, file1: Path, file2: Path) -> tuple[list[TeloType], list[TeloType]]:
    f1_arms: list[TeloType]
    f2_arms: list[TeloType]

    def filter_row(r: dict) -> bool:
        """Returns true if a row should be included."""
        return not (r.get("tvr_len") == "0")  # Flexible TSV format: if tvr_len col doesn't exist, don't crash.

    with open(file1, "r") as fh1, open(file2, "r") as fh2:
        r1 = csv.DictReader(fh1, delimiter="\t")
        r2 = csv.DictReader(fh2, delimiter="\t")

        row: dict
        f1_arms = [scoring.build_telo_from_row(row) for row in filter(filter_row, r1)]
        f2_arms = [scoring.build_telo_from_row(row) for row in filter(filter_row, r2)]

    return f1_arms, f2_arms


def _compute_matrix(scoring: BaseScoringSystem, f1_arms: list[TeloType], f2_arms: list[TeloType]) -> list[list[float]]:
    matrix: list[list[float]] = [[0.0 for _j in f1_arms] for _i in f2_arms]

    for i, f1a in enumerate(f1_arms):
        for j, f2a in enumerate(f2_arms):
            score, cigar = scoring.score_seqs(f1a["tvr_consensus_encoded"], f2a["tvr_consensus_encoded"])
            matrix[j][i] = score
            if score > scoring.score_log_threshold:
                logger.info(
                    f"Found score >{scoring.score_log_threshold}: {_fmt_allele(f1a)} against {_fmt_allele(f2a)}; "
                    f"score: {score:.3f}"
                )
                logger.info(
                    f"  Alignment: \n"
                    f"{_fmt_alignment(f1a['tvr_consensus_encoded'], f2a['tvr_consensus_encoded'], cigar)}"
                )

    return matrix


def _write_outfile(
    f1_arms: list[TeloType],
    f2_arms: list[TeloType],
    out_file: Path | StringIO,
    matrix: list[list[float]],
):
    def _write(f: TextIO):
        header = "\t".join(["", *(_fmt_allele(f1a) for f1a in f1_arms)])
        f.write(f"{header}\n")
        for j, f2a in enumerate(f2_arms):
            f.write("\t".join([_fmt_allele(f2a), *map(str, matrix[j])]) + "\n")

    if isinstance(out_file, StringIO):
        _write(out_file)
    else:
        with open(out_file, "w") as fh:
            _write(fh)


def compare_samples(scoring: BaseScoringSystem, file1: Path, file2: Path, out_file: Path | StringIO):
    # Step 1: load telomere arms from Telogator2/similar TSV files
    f1_arms, f2_arms = _read_sample_files(scoring, file1, file2)
    # Step 2: calculate all-all scoring matrix using TVR sequences
    matrix: list[list[float]] = _compute_matrix(scoring, f1_arms, f2_arms)
    # Step 3: write scoring matrix to out_file
    _write_outfile(f1_arms, f2_arms, out_file, matrix)
