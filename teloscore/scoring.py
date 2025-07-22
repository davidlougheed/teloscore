import math
import parasail

from abc import ABC, abstractmethod

from .types import TeloType
from .utils import decode_cigar

__all__ = [
    "BaseScoringSystem",
    "Scoring1",
    "Scoring2",
]


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
