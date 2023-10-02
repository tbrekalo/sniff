import argparse
import pathlib
import sys

import joblib
import lightgbm as lgb
import polars as pl

FEATURE_COL_NAMES = [
    'query_score',
    'query_lhs_overhang',
    'query_rhs_overhang',

    'target_score',
    'target_lhs_overhang',
    'target_rhs_overhang',
]


def read_sniff_csv(path: pathlib.Path | str) -> pl.DataFrame:
    """Loads sniff output in CSV format into polars data frame.

    Args:
        path: A path towards sniff output

    Returns:
        Polars DataFrame object with loaded data and aditional columns.
    """

    def make_ovlp_expressions(seq_name):
        return [
            (
                (pl.col(f'{seq_name}_end') - pl.col(f'{seq_name}_start')) /
                pl.col(f'{seq_name}_length')
            ).alias(f'{seq_name}_score'),
            (
                pl.col(f'{seq_name}_start') / pl.col(f'{seq_name}_length')
            ).alias(f'{seq_name}_lhs_overhang'),
            (
                (pl.col(f'{seq_name}_length') - pl.col(f'{seq_name}_end')) /
                pl.col(f'{seq_name}_length')
            ).alias(f'{seq_name}_rhs_overhang')
        ]

    return pl.read_csv(
        path
    ).with_columns(
        *make_ovlp_expressions('query'),
        *make_ovlp_expressions('target'),
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='lgbm_filter',
        description='process sniff overlap output and filter correct overlaps',
    )

    parser.add_argument(
        '-m', '--model', type=str, required=True,
        help='path to lgbm model file',
    )

    parser.add_argument(
        '-o', '--overlaps', type=str, required=True,
        help='path to sniff overlaps output in csv format',
    )

    parser.add_argument(
        '-s', '--pair-separator', type=str, default=',',
        help='output pair separator string',
    )

    args = parser.parse_args()

    model = joblib.load(args.model)
    df_overlaps = read_sniff_csv(args.overlaps)

    df_overlaps.with_columns(
        pl.Series(
            model.predict(
                df_overlaps.select(FEATURE_COL_NAMES)
            )
        ).alias('score')
    ).select([
        'query_name', 'target_name', 'score'
    ]).groupby(
        'query_name'
    ).agg(
        pl.all().sort_by('score').last()
    ).groupby(
        'target_name'
    ).agg(
        pl.all().sort_by('score').last()
    ).filter(
        pl.col('score') == 1
    ).select(
        pl.min(pl.col('query_name'), pl.col('target_name')).alias('lhs'),
        pl.max(pl.col('query_name'), pl.col('target_name')).alias('rhs'),
    ).write_csv(
        file=sys.stdout,
        has_header=False,
        separator=args.pair_separator,
    )
