# TeloScore

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14984162.svg)](https://doi.org/10.5281/zenodo.14984162)
[![PyPI version](https://badge.fury.io/py/teloscore.svg)](https://badge.fury.io/py/teloscore)

A telomeric allele comparison-scoring tool for output from [Telogator2](https://github.com/zstephens/telogator2).

Work done on this project is part of a preprint by 
[Zhou *et al.*](https://www.biorxiv.org/content/10.1101/2025.10.07.680721v1.abstract).
If you use this work, please cite this preprint.


## Requirements

* Python 3.10+
* Operating system: tested on macOS and Linux


## Installation

TeloScore can be installed quickly from PyPI using the following command:

```bash
pip install teloscore
```

TeloScore's dependencies are listed in [`pyproject.toml`](./pyproject.toml).


## Usage

The following is an example showing how TeloScore is to be used. The comparison
function should finish in under a minute.

```bash
teloscore compare ./sample_01_child.tsv ./sample_02_parent_1.tsv ./out.tsv
```

> **Note:** There are currently two scoring systems: `1` and `2`, corresponding
> to the scoring systems found in the v0.1.x release and the v0.2.x release.
> The default is `2`. You may specify which scoring system you wish to use
> by passing `--scoring #`, where `#` is `1` or `2`.

The output matrix can then be plotted using the following command:

```bash
teloscore plot ./out.tsv
```


## Running with Example Data

This repository includes Telogator2 telomere calls for the Genome-in-a-Bottle
Ashkenazi trio (HG002-4) in the [`test_data`](./test_data) directory, generated
using reads obtained from https://downloads.pacbcloud.com/public/revio/2022Q4/.

To run the two parental comparisons with these data using scoring system "1"
(the scoring system used in the 
[Zhou *et al.* preprint](https://www.biorxiv.org/content/10.1101/2025.10.07.680721v1.abstract)),
run the following commands:

```bash
# paternal
teloscore compare --scoring 1 \
  ./test_data/tlens_by_allele_HG002.tsv \
  ./test_data/tlens_by_allele_HG003.tsv \
  pat.out.tsv
  
# maternal
teloscore compare --scoring 1 \
  ./test_data/tlens_by_allele_HG002.tsv \
  ./test_data/tlens_by_allele_HG004.tsv \
  mat.out.tsv

# You will now have two files, pat.out.tsv and mat.out.tsv. These can now be 
# plotted using the `teloscore plot` subcommand:
teloscore plot pat.out.tsv
teloscore plot mat.out.tsv
```


## Input Format

TeloScore's `compare` subcommand is primarily designed to accept the 
[TSV output file from Telogator2](https://github.com/zstephens/telogator2?tab=readme-ov-file#output-files).
However, any TSV file with the following columns would work:

* `#chr` / `chr`: Chromosome arm
* `allele_id`: ID for this specific allele (see Telogator2)
* `tvr_consensus`: Allele TVR consensus sequence, in the Telogator2 base encoding 
  (https://github.com/zstephens/telogator2/blob/main/resources/kmers.tsv)


## Output Format

**Definition:** TVR = telomere variable region

The output is a TSV matrix with the TVR of the first sample (child) across the 
columns, and those of the second sample (parent) across the rows. Each entry in
the matrix is a floating-point number between 0 and 1, representing a custom 
similarity score between the two TVRs.


## Copyright Notice

TeloScore is a telomeric allele comparison-scoring tool for output from 
[Telogator2](https://github.com/zstephens/telogator2).

Copyright (C) 2024-2025  McGill University, David Lougheed

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
