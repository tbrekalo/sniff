import argparse
import pathlib

import polars as pl

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


def load_rc_pairs_as_df(pairs_csv: str | pathlib.Path) -> pl.DataFrame:
    try:
        return pl.read_csv(pairs_csv, has_header=False, new_columns=PAIR_COLS)
    except Exception:
        return pl.DataFrame(schema={'query_name': str, 'target_name': str})


def load_paf_ref_mapping_as_df(
        ref_mapping: str | pathlib.Path,
        read_min_coverage: float = 0.875) -> pl.DataFrame:
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
        df_rc_pairs: pl.DataFrame,
        df_ref_reads_mapping: pl.DataFrame) -> pl.DataFrame:
    return df_rc_pairs.join(
        df_ref_reads_mapping.select(
            pl.col('seq_name').alias('query_name'),
            pl.col('length').alias('query_length'),
            pl.col('strand').alias('query_strand'),
            pl.col('ref_start').alias('query_ref_start'),
            pl.col('ref_end').alias('query_ref_end')
        ),
        on='query_name', how='left',
    ).join(
        df_ref_reads_mapping.select(
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='rc_stats',
        description='expand reverse complement pairs with '
        'reference mapping info'
    )

    parser.add_argument(
        '-m', '--mapping', type=str, required=True,
        help='read reference mapping paf',
    )

    parser.add_argument(
        '-p', '--rc-pairs', type=str, required=True,
        help='path to reverse complement csv file',
    )

    args = parser.parse_args()
    df_rc_pairs = load_rc_pairs_as_df(args.rc_pairs)
    df_ref_reads_mapping = load_paf_ref_mapping_as_df(ref_mapping=args.mapping)

    print(expand_pairs_with_mapping(
        df_rc_pairs=df_rc_pairs,
        df_ref_reads_mapping=df_ref_reads_mapping,
    ).write_csv())
