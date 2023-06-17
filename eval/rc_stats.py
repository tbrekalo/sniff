import argparse
import pathlib
from time import perf_counter
from typing import List, Tuple

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

PAF_COLUMNS = [
    'query_name',
    'query_length',
    'query_start',
    'query_end',
    'strand',
    'target_name',
    'target_length',
    'target_start',
    'target_end',
    'n_residue_matches',
]

PAIR_COLS = [
    'query_name',
    'target_name',
]


def load_pairs_df(pairs_csv: str) -> pl.DataFrame:
    try:
        return pl.read_csv(pairs_csv, has_header=False, new_columns=PAIR_COLS)
    except Exception:
        return pl.DataFrame(schema={'query_name': str, 'target_name': str})


def load_paf_ref_mapping(ref_mapping: str, read_min_coverage: float = 0.875) -> pl.DataFrame:
    return pl.read_csv(ref_mapping, separator='\t',
                       has_header=False, columns=list(range(len(PAF_COLUMNS))),
                       new_columns=PAF_COLUMNS).select([
                           pl.col('^query_.*$'),
                           pl.col('strand'),
                           pl.col('target_start').alias('ref_start'),
                           pl.col('target_end').alias('ref_end'),
                       ]).groupby(
        'query_name'
    ).agg(
        pl.all().sort_by(
            pl.col('query_end')-pl.col('query_start') / pl.col('query_length')
        ).last()
    ).filter(
        (pl.col('query_end') - pl.col('query_start')
         ) > pl.col('query_length') * read_min_coverage
    ).select(
        pl.col('query_name').alias('seq_name'),
        pl.col('query_length').alias('length'),
        pl.col('query_start').alias('start'),
        pl.col('query_end').alias('end'),
        pl.col('strand'), pl.col('ref_start'), pl.col('ref_end')
    )


def expand_pairs_with_mapping(
        pairs_df: pl.DataFrame,
        mapping_df: pl.DataFrame) -> pl.DataFrame:
    return pairs_df.join(
        mapping_df.select(
            pl.col('seq_name').alias('query_name'),
            pl.col('length').alias('query_length'),
            pl.col('strand').alias('query_strand'),
            pl.col('ref_start').alias('query_ref_start'),
            pl.col('ref_end').alias('query_ref_end')
        ),
        on='query_name', how='left',
    ).join(
        mapping_df.select(
            pl.col('seq_name').alias('target_name'),
            pl.col('length').alias('target_length'),
            pl.col('strand').alias('target_strand'),
            pl.col('ref_start').alias('target_ref_start'),
            pl.col('ref_end').alias('target_ref_end')
        ),
        on='target_name', how='left',
    ).with_columns(pl.when(
        pl.col('query_strand').is_not_null() &
        pl.col('target_strand').is_not_null()
    ).then(
        pl.max(
            pl.min(pl.col('query_ref_end'), pl.col('target_ref_end')) -
            pl.max(pl.col('query_ref_start'), pl.col('target_ref_start')),
            0
        ) / (
            pl.max(pl.col('query_ref_end'), pl.col('target_ref_end')) -
            pl.min(pl.col('query_ref_start'), pl.col('target_ref_end'))
        )).otherwise(
            0
    ).alias('reference_overlap'),
        pl.when(
        pl.col('query_strand') == pl.col('target_strand')
    ).then('+').otherwise(
        '-'
    ).alias('strand')
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
        ref_mapping_df: pl.DataFrame) -> Tuple[pl.DataFrame, pl.DataFrame]:
    PAIRS_CSV_PATH = '/tmp/sniff-pairs.csv'
    with open(PAIRS_CSV_PATH, 'w+', encoding='UTF-8') as pairs_csv:
        with Popen(create_sniff_spawn_list(
                sniff_path, sniff_args, reads_path), stdout=pairs_csv) as proc:
            peak_memory = 0
            time_start = time_end = perf_counter()

            while proc.poll() is None:
                curr_mem = proc.memory_info().rss
                time_end = perf_counter()

                if curr_mem is not None and curr_mem > peak_memory:
                    peak_memory = curr_mem

    return (
        pl.DataFrame(sniff_args.dict() | RunInfo(
            runtime_s=int(time_end-time_start),
            peak_memory_gib=peak_memory / (2 ** 30)).dict()),
        expand_pairs_with_mapping(
            pairs_df=load_pairs_df(PAIRS_CSV_PATH),
            mapping_df=ref_mapping_df,
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='rc_stats',
        description='calculate precision recall and accuracy for reverse complement pairs'
    )

    parser.add_argument(
        '-m', '--mapping', type=str, required=True,
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
        '-o', '--out', type=str, default='eval-sniff',
        help='path to output folder',
    )

    args = parser.parse_args()
    ref_mapping_df = load_paf_ref_mapping(ref_mapping=args.mapping)
    run_df, mapped_pairs = run_sniff(
        sniff_path=args.sniff_path,
        sniff_args=DEFAULT_ARGS,
        reads_path=args.reads_path,
        ref_mapping_df=ref_mapping_df,
    )

    out_dir = pathlib.Path(args.out)
    if not out_dir.exists():
        out_dir.mkdir()

    out_run_info = out_dir.joinpath('run-info.csv')
    out_pairs_info = out_dir.joinpath('pairs-mapped.csv')

    run_df.write_csv(out_run_info)
    mapped_pairs.write_csv(out_pairs_info)
