from typing import Iterable

__all__ = ["decode_cigar"]

# Taken from STRkit: https://github.com/davidlougheed/strkit/blob/master/strkit/call/cigar.py


def _decode_cigar_item(item: int) -> tuple[int, int]:
    return item & 15, item >> 4


def decode_cigar(encoded_cigar: list[int]) -> Iterable[tuple[int, int]]:
    return map(_decode_cigar_item, encoded_cigar)


# End taken from STRkit
