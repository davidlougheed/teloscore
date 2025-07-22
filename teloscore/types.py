from typing import TypedDict

__all__ = ["TeloType"]


class TeloType(TypedDict):
    arm: str
    allele_id: str  # used to be an int, but we can also have ^\d\di$ pattern allele IDs
    tvr_consensus: str
    tvr_consensus_encoded: str
