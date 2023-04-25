# Sniff

Heuristic tool for pairing up reverse complement reads in fasta/fastq files.

## Build dependencies

- gcc 9.4 or higher
- clang 10.1 or higher
- intel tbb 2020.3
  - should be compatible with later oneapi versions
- cmake 3.21 or higher
- git 2.25.1 or higher
  - earlier version should do just fine
  - git is required for cmake to fetch par of internal dependencies

## Usage

```bash
pair up potential reverse complement reads
Usage:
  sniff [OPTION...] positional parameters

 general options:
  -h, --help         print help
  -v, --version      print version
  -t, --threads arg  number of threads to use (default: 1)

 heuristic options:
  -p, --percent arg        maximum allowed difference in length as % of
                           shorter read's length (default: 0.05)
  -l, --sample-length arg  maximum sample length from beginning/end of
                           sequence (default: 5000)
  -e, --edit-distance arg  maximum allowed edit distance between samples
                           (default: 100)

 input options:
      --input arg  input fasta/fastq file
```

## Methods

For each read `sniff` finds reads which are similar in length and calculates edit-distance score between first K bases of target read and K bases of reverse complement from the similar in length read. If the computed edit-distance is below a certain threshold, reads are marked as a reverse complement pair.

## Evaluation

TODO...
