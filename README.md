# Sniff

Heuristic tool for pairing up reverse complement reads in fasta/fastq files.

## Methods

Sniff loads sequences in batches displaying progress along the way. Once all sequences have been loaded into memory they are again processed in batches. Each batch is used for construcing a target index from reverse complemented reads. Original reads are then mapped against the constructed index
using a seed and chain approach with minor modifications. For each read we remember the strongest matching reverse complement read and output reverse complement pairs. Each read appears in at most one pair.

## Usage

After building:

```bash
pair up potential reverse complement reads
Usage:
  sniff [OPTION...] <reads>

 general options:
  -h, --help         print help
  -v, --version      print version
  -t, --threads arg  number of threads to use (default: 1)

 heuristic options:
  -a, --alpha arg  shorter read length as percentage of longer read lenght 
                   in pair (default: 0.10)
  -b, --beta arg   minimum required coverage on each read (default: 0.98)

 input options:
      --input arg  input fasta/fastq file

 mapping options:
  -k, --kmer-length arg    kmer length used in mapping (default: 15)
  -w, --window-length arg  window length used in mapping (default: 5)
  -f, --frequent arg       filter f most frequent kmers (default: 0.0002)

```

## Dependencies

- linux kernel 2.6.32 or higher
- gcc 10.4 or higher
- clang 10.1 or higher
- intel tbb 2020.3
  - should be compatible with later oneapi versions
- cmake 3.21 or higher
- git 2.25.1 or higher
  - earlier version should do just fine
  - git is required for cmake to fetch par of internal dependencies

### Test (optional) dependencies

- Catch2
  - fetched via cmake if missing

## Build

After running running git clone:

```bash
git clone git@github.com:tbrekalo/sniff.git
```

Run the following:

```bash
cmake -S ./sniff -B ./sniff/build -DCMAKE_BUILD_TYPE=Release -G Ninja
cmake --build build
```

You can ommit the `-G Ninja` option in case your host system doesn't support the `Ninja` build system.
