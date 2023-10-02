# Sniff

Heuristic tool for pairing up reverse complement reads in fasta/fastq files.

![ci](https://github.com/tbrekalo/sniff/actions/workflows/ci.yml/badge.svg?branch=master)

## Methods

Sniff loads sequences in batches displaying progress along the way. Once all sequences have been loaded into memory they are again processed in batches. Each batch is used for construcing a target index from reverse complemented reads. Original reads are then mapped against the constructed index using a seed and chain approach with minor modifications. For each read we remember the strongest matching reverse complement read and output reverse complement pairs as overlaps. Here we define an overlap as a tuple `query_name, query_start, query_end, target_name, target_start, target_end`. Later those overlaps are processed with a pretrained machine learning model outputing the final result in a csv/tsv format for later use.

## Dependencies

### C++
- linux kernel 2.6.32 or higher
- gcc 11 or higher
- clang 11 or higher
- intel tbb 2020.3
  - should be compatible with later oneapi versions
- cmake 3.21 or higher
- git 2.25.1 or higher
  - earlier version should do just fine
  - git is required for cmake to fetch par of internal dependencies

#### Test (optional) dependencies

- Catch2
  - fetched via cmake if missing

### Python
conan2==0.0.4
joblib==1.3.2
lightgbm==4.1.0
polars==0.18.5
psutil==5.9.5
pydantic==1.10.9
scikit-learn==1.3.1

## Build

```bash
git clone git@github.com:tbrekalo/sniff.git
cd sniff
make release
```

## Usage

From sniff root directory:

```bash
source ./venv/bin/activate
sniff -t 32 path_to_reads.fasta > /tmp/sniff-ovlps.csv
python ./scripts/inference/lgbm_filter.py -m resources/sniff-lgbm-model.pkl -o /tmp/sniff.csv > pairs.csv
```
