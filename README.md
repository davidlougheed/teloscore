# `comptg2`

A telomeric allele comparison-scoring tool for output from [Telogator2](https://github.com/zstephens/telogator2).


## Installation

`comptg2` can be installed from PyPI using the following command:

```bash
pip install comptg2
```


## Usage

The following is an example showing how `comptg2` is to be used:

```bash
comptg2 compare ./sample_01_child.tsv ./sample_02_parent_1.tsv ./out.tsv
```

The output matrix can then be plotted using the following command:

```bash
comptg2 plot ./out.tsv
```


## Output Format

**Definition:** TVR = telomere variable region

The output is a TSV matrix with the TVR of the first sample (child) across the 
columns, and those of the second sample (parent) across the rows. Each entry in
the matrix is a floating-point number between 0 and 1, representing a custom 
similarity score between the two TVRs.


## Copyright Notice

`comptg2` is a telomeric allele comparison-scoring tool for output from 
[Telogator2](https://github.com/zstephens/telogator2).

Copyright (C) 2024  McGill University, David Lougheed

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
