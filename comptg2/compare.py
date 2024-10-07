import csv
import logging
import math
import parasail

from typing import Iterable, TypedDict

__all__ = [
    "re_encode_log",
    "compare_samples",
]


GAP_OPEN_PENALTY = 3
GAP_EXTEND_PENALTY = 2
MATCH_SCORE = 5
CANONICAL_MATCH_SCORE = 1
CANONICAL_VS_INDEL_MISMATCH_SCORE = -1  # between negative mismatch score and positive canonical match score...
DUBIOUS_MATCH_SCORE = 2
MISMATCH_SCORE = -4

COMPRESSION_LOG_BASE = 4

# TODO: could maybe implement a more advanced scoring matrix
# Alphabet is from https://github.com/zstephens/telogator2/blob/main/resources/kmers.tsv
# 'A' is 'UNKNOWN_LETTER'
SCORING_ALPHABET = "ACDEFGHIKLMNPRSQTVWY"
CANONICAL_LETTER = "C"
CANONICAL_INDEX = SCORING_ALPHABET.index(CANONICAL_LETTER)
DUBIOUS_LETTERS = "VWY"

SCORING_MATRIX = parasail.matrix_create(SCORING_ALPHABET, MATCH_SCORE, MISMATCH_SCORE)

MATCH_SCORES = {
    tuple(CANONICAL_LETTER): CANONICAL_MATCH_SCORE,
    tuple(DUBIOUS_LETTERS): MATCH_SCORE,
    tuple(set(SCORING_ALPHABET) - {CANONICAL_LETTER} - set(DUBIOUS_LETTERS)): MATCH_SCORE,
}

# There are a lot of canonical motifs
#  - give them a lower match score to weight importance to patterns of non-canonical motifs.
SCORING_MATRIX[CANONICAL_INDEX, CANONICAL_INDEX] = CANONICAL_MATCH_SCORE

CANONICAL_INDEL_INDEX = SCORING_ALPHABET.index("T")
CANONICAL_INDEL_DUBIOUS_INDEX = SCORING_ALPHABET.index("V")

# Give canonical letters matching with canonical indel variant a lesser penalty, since these types of errors should be
# more likely to be sequencing error.
SCORING_MATRIX[CANONICAL_INDEX, CANONICAL_INDEL_INDEX] = CANONICAL_VS_INDEL_MISMATCH_SCORE
SCORING_MATRIX[CANONICAL_INDEL_INDEX, CANONICAL_INDEX] = CANONICAL_VS_INDEL_MISMATCH_SCORE
# TODO: dubious

# Give dubious letter self-matches a lower positive score, to reduce their influence on the final score.
for d in DUBIOUS_LETTERS:
    SCORING_MATRIX[SCORING_ALPHABET.index(d), SCORING_ALPHABET.index(d)] = DUBIOUS_MATCH_SCORE


logging.basicConfig(level=logging.INFO)

logger = logging.getLogger("comptg2")


# Taken from STRkit: https://github.com/davidlougheed/strkit/blob/master/strkit/call/cigar.py


def _decode_cigar_item(item: int) -> tuple[int, int]:
    return item & 15, item >> 4


def decode_cigar(encoded_cigar: list[int]) -> Iterable[tuple[int, int]]:
    return map(_decode_cigar_item, encoded_cigar)


# End taken from STRkit


def match_score_for_letter(x: str):
    for k, v in MATCH_SCORES.items():
        if x in k:
            return v


class TeloType(TypedDict):
    arm: str
    allele_id: str  # used to be an int, but we can also have ^\d\di$ pattern allele IDs
    tvr_consensus: str
    tvr_consensus_encoded: str


def re_encode_log(seq: str) -> str:
    assert len(seq) > 0

    revised_seq: list[str] = []

    current_char: str = seq[0]
    current_count: int = 1

    for c in (*seq[1:], ""):
        if c == current_char:
            current_count += 1
        else:
            revised_seq.append("".join(current_char * int(1.0 + round(math.log(current_count, COMPRESSION_LOG_BASE)))))

            current_char = c
            current_count = 1

    return "".join(revised_seq)


def score_seqs(seq1: str, seq2: str) -> tuple[float, tuple[tuple[int, int], ...]]:
    qs, dbs = (seq1, seq2) if len(seq1) <= len(seq2) else (seq2, seq1)

    r = parasail.nw_trace_striped_32(qs, dbs, GAP_OPEN_PENALTY, GAP_EXTEND_PENALTY, SCORING_MATRIX)
    cigar = r.cigar

    total_possible_score: int = sum(map(match_score_for_letter, qs))
    final_score: float = max(r.score / total_possible_score, 0.0)
    return final_score, tuple(decode_cigar(cigar.seq))


def build_telo_from_row(row: dict) -> TeloType:
    tvr = row["tvr_consensus"]
    return {
        "arm": row["#chr"],
        "allele_id": row["allele_id"],
        "tvr_consensus": tvr,
        "tvr_consensus_encoded": re_encode_log(tvr),
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


def compare_samples(file1: str, file2: str, out_file: str):
    f1_arms: list[TeloType] = []
    f2_arms: list[TeloType] = []

    with open(file1, "r") as fh1, open(file2, "r") as fh2:
        r1 = csv.DictReader(fh1, delimiter="\t")
        r2 = csv.DictReader(fh2, delimiter="\t")

        row: dict
        for row in r1:
            if row["tvr_len"] == "0":
                continue
            f1_arms.append(build_telo_from_row(row))
        for row in r2:
            if row["tvr_len"] == "0":
                continue
            f2_arms.append(build_telo_from_row(row))

    matrix: list[list[float]] = [[0.0 for _j in f1_arms] for _i in f2_arms]

    for i, f1a in enumerate(f1_arms):
        for j, f2a in enumerate(f2_arms):
            score, cigar = score_seqs(f1a["tvr_consensus_encoded"], f2a["tvr_consensus_encoded"])
            matrix[j][i] = score
            if score > 0.6:
                logger.info(f"Found score >0.6: {_fmt_allele(f1a)} against {_fmt_allele(f2a)}; score: {score:.3f}")
                logger.info(
                    f"  Alignment: \n"
                    f"{_fmt_alignment(f1a['tvr_consensus_encoded'], f2a['tvr_consensus_encoded'], cigar)}"
                )

    with open(out_file, "w") as fh:
        header = "\t".join(["", *(_fmt_allele(f1a) for f1a in f1_arms)])
        fh.write(f"{header}\n")
        for j, f2a in enumerate(f2_arms):
            fh.write("\t".join([_fmt_allele(f2a), *map(str, matrix[j])]) + "\n")
