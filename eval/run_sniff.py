import argparse
import pathlib
from time import perf_counter
from typing import List

import polars as pl
from psutil import Popen
from pydantic import BaseModel


class SniffArgs(BaseModel):
    threads: int
    percent: float
    kmer_length: int
    window_length: int
    chain: int
    gap: int


class RunInfo(BaseModel):
    runtime_s: float
    peak_memory_gib: float


DF_COLS = [
    'threads',
    'percent',
    'kmer_length',
    'window_length',
    'chain',
    'gap',
    'runtime_s',
    'peak_memory_gib',
    'recall',
    'precision',
]

DEFAULT_ARGS = SniffArgs(
    threads=32,
    percent=0.01,
    kmer_length=15,
    window_length=5,
    chain=4,
    gap=500,
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
        sniff_path: str, sniff_args: SniffArgs, reads_path: str) -> List[str]:
    return [sniff_path, *format_sniff_args(sniff_args, reads_path)]


def run_sniff(
        sniff_path: str,
        sniff_args: SniffArgs,
        reads_path: str,
        pairs_out_path: str | pathlib.Path) -> pl.DataFrame:
    with open(pairs_out_path, 'w+', encoding='UTF-8') as pairs_csv:
        with Popen(create_sniff_spawn_list(
                sniff_path, sniff_args, reads_path), stdout=pairs_csv) as proc:
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
        '-o', '--out', type=str, default='./eval-sniff',
        help='path to output folder',
    )

    args = parser.parse_args()
    out_dir = pathlib.Path(args.out)
    if not out_dir.exists():
        out_dir.mkdir()

    run_sniff(
        sniff_path=args.sniff_path,
        sniff_args=DEFAULT_ARGS,
        reads_path=args.reads_path,
        pairs_out_path=out_dir.joinpath('rc-pairs.csv'),
    ).write_csv(out_dir.joinpath('run-info.csv'))
