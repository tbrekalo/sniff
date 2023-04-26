import random
from typing import Dict, List, Tuple

import numba
import numpy as np

BASES = 'ACGT'


@numba.njit
def encode_kmer(kmer: str) -> int:
    MAPPING: Dict[str, int] = {
        'A': 0,
        'C': 1,
        'G': 2,
        'T': 3,
    }

    dst = 0
    for b in kmer:
        dst = (dst << 2) | MAPPING[b]
    return dst


@numba.njit
def hash(key):
    key = np.uint64(key)
    # key = (key << 21) - key - 1;
    key = np.uint64((~key) + (key << np.uint64(21)))
    key = np.uint64(key ^ (key >> np.uint64(24)))
    key = np.uint64((key + (key << np.uint64(3))) +
                    (key << np.uint64(8)))  # key * 265
    key = np.uint64(key ^ (key >> np.uint64(14)))
    key = np.uint64((key + (key << np.uint64(2))) +
                    (key << np.uint64(4)))  # key * 21
    key = np.uint64(key ^ (key >> np.uint64(28)))
    key = np.uint64(key + (key << np.uint64(31)))
    return key


@numba.njit
def find_kmers(data: str, k: int) -> List[Tuple[int, np.uint64, str]]:
    return [
        (i, hash(encode_kmer(data[i:i+k])), data[i:i+k])
        for i in range(0, max(0, len(data) - k))
    ]


def pick_minimizers(data: str, k: int, w: int) -> List[Tuple[int, np.uint64, str]]:
    kmers = find_kmers(data, k)
    return sorted(list(
        set((
            min(kmers[i:i+w], key=lambda x: x[1])
            for i in range(0, len(kmers) - w)
            if i + w < len(data)
        ))
    ))


SAMPLE = 'CTTTTCAAATATATG'
for t in find_kmers(SAMPLE, 5):
    print(t)
print()
for m in pick_minimizers(SAMPLE, 5, 5):
    print(m)
