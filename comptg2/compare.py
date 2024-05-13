import csv
import math
import parasail

from typing import TypedDict

__all__ = [
    "re_encode_log",
    "compare_samples",
]

GAP_OPEN_PENALTY = 3
GAP_EXTEND_PENALTY = 1
MATCH_SCORE = 5
MISMATCH_SCORE = -4

COMPRESSION_LOG_BASE = 5

# TODO: could maybe implement a more advanced scoring matrix
# Alphabet is from https://github.com/zstephens/telogator2/blob/main/resources/kmers.tsv
# 'A' is 'UNKNOWN_LETTER'
SCORING_MATRIX = parasail.matrix_create("ACDEFGHIKLMNPRSQTVWY", MATCH_SCORE, MISMATCH_SCORE)


class TeloType(TypedDict):
    arm: str
    allele_id: int
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


def score_seqs(seq1: str, seq2: str) -> float:
    qs, dbs = (seq1, seq2) if len(seq1) <= len(seq2) else (seq2, seq1)

    r = parasail.nw_trace_striped_32(qs, dbs, GAP_OPEN_PENALTY, GAP_EXTEND_PENALTY, SCORING_MATRIX)

    # cigar = r.cigar
    # print(seq1, seq2, cigar.decode, r.score)

    final_score: float = max(r.score / (MATCH_SCORE * len(qs)), 0.0)
    return final_score


def build_telo_from_row(row: dict) -> TeloType:
    tvr = row["tvr_consensus"]
    return {
        "arm": row["#chr"],
        "allele_id": int(row["allele_id"]),
        "tvr_consensus": tvr,
        "tvr_consensus_encoded": re_encode_log(tvr),
    }


def compare_samples(file1: str, file2: str):
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

    for f1a in f1_arms:
        for f2a in f2_arms:
            # print("---------")
            score = score_seqs(f1a["tvr_consensus_encoded"], f2a["tvr_consensus_encoded"])
            if score > 0.8:
                print(
                    "passed",
                    "1",
                    f1a["arm"],
                    f"({f1a['allele_id']})",
                    "2",
                    f2a["arm"],
                    f"({f2a['allele_id']})",
                    "score:",
                    score,
                )
                print(f1a["tvr_consensus_encoded"])
                print(f2a["tvr_consensus_encoded"])
