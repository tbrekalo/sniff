''' script for generating fasta test data '''

import argparse
import datetime
import random

BASES = 'ACGT'
COMPLEMENTS = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        'gen_seq', 'generate DNA sequence of given length')
    parser.add_argument('len', type=int, help='sequence length')
    parser.add_argument(
        '-r', '--reverse-complement',
        action='store_true', help='print reverse complement sequence')

    args = parser.parse_args()

    seq_name = '>' + datetime.datetime.now().strftime(
        '%d-%m-%y-%H-%M-%S'
    )

    seq_data = ''.join(
        (random.choice(BASES) for _ in range(args.len))
    )

    print(f'{seq_name}\n{seq_data}')
    if args.reverse_complement:
        rc_seq_name = seq_name + '_rc'
        rc_seq_data = ''.join(
            COMPLEMENTS[b] for b in seq_data[::-1]
        )

        print(f'{rc_seq_name}\n{rc_seq_data}')
