import argparse
import pathlib
import sys
from time import perf_counter
from typing import List, Tuple

import polars as pl
from psutil import Popen
from pydantic import BaseModel


class SniffArgs(BaseModel):
    threads: int
    alpha_percent: float
    beta_percent: float
    kmer_length: int
    window_length: int
    chain: int
    gap: int
    frequent: float


class RunInfo(BaseModel):
    runtime_s: float
    peak_memory_gib: float


DF_COLS = [
    'threads',
    'alpha_percent',
    'beta_percent',
    'kmer_length',
    'window_length',
    'chain',
    'gap',
    'frequent',
    'runtime_s',
    'peak_memory_gib',
    'recall',
    'precision',
]

DEFAULT_ARGS = SniffArgs(
    threads=32,
    alpha_percent=0.1,
    beta_percent=0.98,
    kmer_length=15,
    window_length=5,
    chain=4,
    gap=500,
    frequent=0.001,
)


def format_sniff_args(sniff_args: SniffArgs, reads_path: str) -> List[str]:
    dst = [
        val for k, v in sniff_args.dict().items()
        for val in ('--' + k.replace('_', '-'), str(v))
        if k != 'minhash'
    ]

    dst.append(reads_path)
    return dst


def create_sniff_spawn_list(
        sniff_path: str | pathlib.Path,
        sniff_args: SniffArgs,
        reads_path: str | pathlib.Path) -> List[str]:
    return [
        str(sniff_path), *format_sniff_args(sniff_args, str(reads_path))
    ]


def run_sniff(
        sniff_path: str,
        sniff_args: SniffArgs,
        reads_path: str | pathlib.Path) -> pl.DataFrame:

    with Popen(create_sniff_spawn_list(
        sniff_path, sniff_args, reads_path)
    ) as proc:
        peak_memory = 0
        time_start = time_end = perf_counter()

        while proc.poll() is None:
            curr_mem = proc.memory_info().rss
            time_end = perf_counter()

            if curr_mem is not None and curr_mem > peak_memory:
                peak_memory = curr_mem

    return pl.DataFrame(sniff_args.dict() | RunInfo(
        runtime_s=int(time_end-time_start),
        peak_memory_gib=peak_memory / (2 ** 30)).dict()
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='rc_stats',
        description='run sniff and record runtime information'
    )

    parser.add_argument(
        '-s', '--sniff-path', type=str, required=True,
        help='sniff path',
    )

    parser.add_argument(
        '-r', '--reads-path', type=str, required=True,
        help='path to reads',
    )

    parser.add_argument(
        '-o', '--out', type=str, nargs='?', default=None,
        help='file to output runtime info; otherwise prints to stderr',
    )

    args = parser.parse_args()
    df_run_info = run_sniff(
        sniff_path=args.sniff_path,
        sniff_args=DEFAULT_ARGS,
        reads_path=args.reads_path,
    )

    if args.out is not None:
        with open(args.out, 'w+', encoding='utf-8') as f:
            print(df_run_info.write_csv(), file=f)
    else:
        print(df_run_info, file=sys.stderr)
