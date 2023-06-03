import argparse
from time import perf_counter
from typing import List, Tuple

import polars as pl
from psutil import Popen
from pydantic import BaseModel


class SniffArgs(BaseModel):
    threads: int
    percent: float
    query_length: int
    kmer_length: int
    window_length: int
    chain: int
    gap: int


class RunInfo(BaseModel):
    runtime_s: float
    peak_memory_gib: float
    recall: float
    precision: float


DF_COLS = [
    'threads',
    'percent',
    'query_length',
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
    query_length=5000,
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


def calc_precision_recall(
        gt_df: pl.DataFrame, paris_csv: str) -> Tuple[float, float]:

    try:
        predicted_df = pl.read_csv(
            paris_csv,
            has_header=False,
            new_columns=['query_name', 'target_name']
        )
    except:
        predicted_df = pl.DataFrame(
            None, [('query_name', str), ('target_name', str)])

    # true positives dataframe
    tp_df = predicted_df.join(
        gt_df,
        on=['query_name', 'target_name']
    )

    return (
        (tp_df.height / predicted_df.height
         if not predicted_df.is_empty() else 0),
        tp_df.height / gt_df.height
    )


def run_sniff(
        sniff_path: str,
        gt_df: pl.DataFrame,
        sniff_args: SniffArgs,
        reads_path: str) -> RunInfo:
    pairs_csv_path = '/tmp/sniff-pairs.csv'
    with open(pairs_csv_path, 'w+', encoding='UTF-8') as pairs_csv:
        with Popen(create_sniff_spawn_list(
                sniff_path, sniff_args, reads_path), stdout=pairs_csv) as proc:
            peak_memory = 0
            time_start = time_end = perf_counter()

            while proc.poll() is None:
                curr_mem = proc.memory_info().rss
                time_end = perf_counter()

                if curr_mem is not None and curr_mem > peak_memory:
                    peak_memory = curr_mem

    precision, recall = calc_precision_recall(gt_df, pairs_csv_path)
    return RunInfo(
        runtime_s=int(time_end-time_start),
        peak_memory_gib=peak_memory / (2 ** 30),
        precision=precision,
        recall=recall,
    )


def grid_serach(sniff_path: str, reads_path: str, gt_df: pl.DataFrame) -> pl.DataFrame:
    runs = []

    sniff_args = DEFAULT_ARGS.copy()
    runs.append(
        sniff_args.dict() |
        run_sniff(sniff_path, gt_df, sniff_args, reads_path).dict()
    )

    return pl.DataFrame(runs)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='rc_stats',
        description='calculate precision recall and accuracy for reverse complement pairs'
    )

    parser.add_argument(
        '-g', '--ground-truth', type=str, required=True,
        help='read reference mapping paf',
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
        '-o', '--out', type=str, default='sniff-eval.csv',
        help='path to output csv',
    )

    args = parser.parse_args()
    ground_truth_df = pl.read_csv(
        args.ground_truth
    ).select(
        pl.min(pl.col('query_name'), pl.col(
            'target_name')).alias('query_name'),
        pl.max(pl.col('query_name'), pl.col(
            'target_name')).alias('target_name'),
    ).filter(
        pl.col('query_name') != pl.col('target_name')
    )

    res_df = grid_serach(
        sniff_path=args.sniff_path,
        reads_path=args.reads_path,
        gt_df=ground_truth_df,
    )

    res_df.write_csv(args.out)
