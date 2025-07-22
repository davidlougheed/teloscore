import csv
import logging
import math
import parasail

from abc import ABC, abstractmethod
from typing import Iterable, TypedDict

__all__ = [
    "BaseScoringSystem",
    "Scoring1",
    "Scoring2",
    "compare_samples",
]


class TeloType(TypedDict):
    arm: str
    allele_id: str  # used to be an int, but we can also have ^\d\di$ pattern allele IDs
    tvr_consensus: str
    tvr_consensus_encoded: str


class BaseScoringSystem(ABC):
    def __init__(self):
        self.SCORING_MATRIX = self.build_scoring_matrix()

    COMPRESSION_LOG_BASE = 4

    SCORE_LOG_THRESHOLD = 0.8

    @abstractmethod
    def build_scoring_matrix(self) -> parasail.Matrix:
        pass

    def encode_seq(self, seq: str) -> str:
        assert len(seq) > 0

        revised_seq: list[str] = []

        current_char: str = seq[0]
        current_count: int = 1

        for c in (*seq[1:], ""):
            if c == current_char:
                current_count += 1
            else:
                revised_seq.append(
                    "".join(current_char * int(1.0 + round(math.log(current_count, self.COMPRESSION_LOG_BASE))))
                )

                current_char = c
                current_count = 1

        return "".join(revised_seq)

    @abstractmethod
    def score_seqs(self, seq1: str, seq2: str) -> tuple[float, tuple[tuple[int, int], ...]]:
        pass

    def build_telo_from_row(self, row: dict) -> TeloType:
        # required columns in row: #chr/chr    allele_id    tvr_consensus
        tvr = row["tvr_consensus"]
        arm = row.get("#chr", row.get("chr", ""))
        if not arm:
            raise ValueError(f"Could not get chromosome arm from row: {row}")
        return {
            "arm": arm,
            "allele_id": row["allele_id"],
            "tvr_consensus": tvr,
            "tvr_consensus_encoded": self.encode_seq(tvr),
        }


class Scoring1(BaseScoringSystem):
    GAP_OPEN_PENALTY = 3
    GAP_EXTEND_PENALTY = 1
    MATCH_SCORE = 5
    CANONICAL_MATCH_SCORE = 1
    MISMATCH_SCORE = -4

    COMPRESSION_LOG_BASE = 5

    # TODO: could maybe implement a more advanced scoring matrix
    # Alphabet is from https://github.com/zstephens/telogator2/blob/main/resources/kmers.tsv
    # 'A' is 'UNKNOWN_LETTER'
    SCORING_ALPHABET = "ACDEFGHIKLMNPRSQTVWY"
    CANONICAL_LETTER = "C"

    SCORE_LOG_THRESHOLD = 0.8

    def build_scoring_matrix(self) -> parasail.Matrix:
        canonical_index = self.SCORING_ALPHABET.index(self.CANONICAL_LETTER)
        m = parasail.matrix_create("ACDEFGHIKLMNPRSQTVWY", self.MATCH_SCORE, self.MISMATCH_SCORE)
        m[canonical_index, canonical_index] = self.CANONICAL_MATCH_SCORE
        return m

    def score_seqs(self, seq1: str, seq2: str) -> tuple[float, tuple[tuple[int, int], ...]]:
        qs, dbs = (seq1, seq2) if len(seq1) <= len(seq2) else (seq2, seq1)

        r = parasail.nw_trace_striped_32(qs, dbs, self.GAP_OPEN_PENALTY, self.GAP_EXTEND_PENALTY, self.SCORING_MATRIX)
        cigar = r.cigar

        total_possible_score: int = sum(
            (self.CANONICAL_MATCH_SCORE if qc == self.CANONICAL_LETTER else self.MATCH_SCORE) for qc in qs
        )
        final_score: float = max(r.score / total_possible_score, 0.0)
        return final_score, tuple(decode_cigar(cigar.seq))


class Scoring2(BaseScoringSystem):
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
    DUBIOUS_LETTERS = "VWY"

    MATCH_SCORES = {
        tuple(CANONICAL_LETTER): CANONICAL_MATCH_SCORE,
        tuple(DUBIOUS_LETTERS): MATCH_SCORE,
        tuple(set(SCORING_ALPHABET) - {CANONICAL_LETTER} - set(DUBIOUS_LETTERS)): MATCH_SCORE,
    }

    SCORE_LOG_THRESHOLD = 0.6

    def build_scoring_matrix(self) -> parasail.Matrix:
        canonical_index = self.SCORING_ALPHABET.index(self.CANONICAL_LETTER)
        canonical_indel_index = self.SCORING_ALPHABET.index("T")

        m = parasail.matrix_create(self.SCORING_ALPHABET, self.MATCH_SCORE, self.MISMATCH_SCORE)

        # There are a lot of canonical motifs
        #  - give them a lower match score to weight importance to patterns of non-canonical motifs.
        m[canonical_index, canonical_index] = self.CANONICAL_MATCH_SCORE

        # Give canonical letters matching with canonical indel variant a lesser penalty, since these types of errors
        # should be more likely to be sequencing error.
        m[canonical_index, canonical_indel_index] = self.CANONICAL_VS_INDEL_MISMATCH_SCORE
        m[canonical_indel_index, canonical_index] = self.CANONICAL_VS_INDEL_MISMATCH_SCORE
        # TODO: dubious?

        # Give dubious letter self-matches a lower positive score, to reduce their influence on the final score.
        for d in self.DUBIOUS_LETTERS:
            m[self.SCORING_ALPHABET.index(d), self.SCORING_ALPHABET.index(d)] = self.DUBIOUS_MATCH_SCORE

        return m

    def match_score_for_letter(self, x: str):
        for k, v in self.MATCH_SCORES.items():
            if x in k:
                return v
        raise ValueError(f"Invalid letter encountered: {x}")

    def score_seqs(self, seq1: str, seq2: str) -> tuple[float, tuple[tuple[int, int], ...]]:
        qs, dbs = (seq1, seq2) if len(seq1) <= len(seq2) else (seq2, seq1)

        r = parasail.nw_trace_striped_32(qs, dbs, self.GAP_OPEN_PENALTY, self.GAP_EXTEND_PENALTY, self.SCORING_MATRIX)
        cigar = r.cigar

        total_possible_score: int = sum(map(self.match_score_for_letter, qs))
        final_score: float = max(r.score / total_possible_score, 0.0)
        return final_score, tuple(decode_cigar(cigar.seq))


logging.basicConfig(level=logging.INFO)

logger = logging.getLogger("teloscore")


# Taken from STRkit: https://github.com/davidlougheed/strkit/blob/master/strkit/call/cigar.py


def _decode_cigar_item(item: int) -> tuple[int, int]:
    return item & 15, item >> 4


def decode_cigar(encoded_cigar: list[int]) -> Iterable[tuple[int, int]]:
    return map(_decode_cigar_item, encoded_cigar)


# End taken from STRkit


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


def _read_sample_files(scoring: BaseScoringSystem, file1: str, file2: str) -> tuple[list[TeloType], list[TeloType]]:
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
            if score > scoring.SCORE_LOG_THRESHOLD:
                logger.info(
                    f"Found score >{scoring.SCORE_LOG_THRESHOLD}: {_fmt_allele(f1a)} against {_fmt_allele(f2a)}; "
                    f"score: {score:.3f}"
                )
                logger.info(
                    f"  Alignment: \n"
                    f"{_fmt_alignment(f1a['tvr_consensus_encoded'], f2a['tvr_consensus_encoded'], cigar)}"
                )

    return matrix


def _write_outfile(f1_arms: list[TeloType], f2_arms: list[TeloType], out_file: str, matrix: list[list[float]]):
    with open(out_file, "w") as fh:
        header = "\t".join(["", *(_fmt_allele(f1a) for f1a in f1_arms)])
        fh.write(f"{header}\n")
        for j, f2a in enumerate(f2_arms):
            fh.write("\t".join([_fmt_allele(f2a), *map(str, matrix[j])]) + "\n")


def compare_samples(scoring: BaseScoringSystem, file1: str, file2: str, out_file: str):
    # Step 1: load telomere arms from Telogator2/similar TSV files
    f1_arms, f2_arms = _read_sample_files(scoring, file1, file2)
    # Step 2: calculate all-all scoring matrix using TVR sequences
    matrix: list[list[float]] = _compute_matrix(scoring, f1_arms, f2_arms)
    # Step 3: write scoring matrix to out_file
    _write_outfile(f1_arms, f2_arms, out_file, matrix)
