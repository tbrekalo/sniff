import argparse
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
def hash(key: np.uint64, kMask: np.uint64):
    key = ((~key) + (key << np.uint64(21))) & kMask
    key = key ^ (key >> np.uint64(24))
    key = ((key + (key << np.uint64(3))) + (key << np.uint64(8))) & kMask
    key = key ^ (key >> np.uint64(14))
    key = ((key + (key << np.uint64(2))) + (key << np.uint64(4))) & kMask
    key = key ^ (key >> np.uint64(28))
    key = (key + (key << np.uint64(31))) & kMask
    return key


@numba.njit
def find_kmers(data: str, k: int) -> List[Tuple[int, np.uint64, str]]:
    mask = np.uint64((np.uint64(1) << (2 * k)) - 1)
    return [
        (i, hash(encode_kmer(data[i:i+k]), mask), data[i:i+k])
        for i in range(0, max(0, len(data) - k))
    ]


def find_minimizers(data: str, k: int, w: int) -> List[Tuple[int, np.uint64, str]]:
    kmers = find_kmers(data, k)
    return sorted(list(
        set((
            min(kmers[i: i+w], key=lambda x: x[1])
            for i in range(0, len(kmers) - w)
            if i + w < len(data)
        ))
    ))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        'minimize', 'prints a set of minimizers for the input string')
    parser.add_argument('sequence', type=str, help='input sequence')
    parser.add_argument('-k', '--kmer-length',
                        type=int, default=15, help='kmer length')
    parser.add_argument('-w', '--window-length',
                        type=int, default=5, help='window length')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='output in tsv with included hash values')

    args = parser.parse_args()
    for pos, hash, mini in find_minimizers(
            args.sequence, args.kmer_length, args.window_length):
        if not args.verbose:
            print(f'{pos},{mini}')
        else:
            print(f'{pos}\t{hash}\t{mini}')
