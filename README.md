# nei_vcf

Calculate variant-based Nei distances using iPyrad's VCF output.

## Installation

```
git clone https://github.com/TankredO/nei_vcf
cd nei_vcf
pip install .
```

## Usage

```
usage: nei_vcf [-o OUT_FILE] [-n N_COL_SKIP] [-e ENTRY_IDX] [-v] [-h] vcf

Calculate variant-based Nei distances using a VCF input file.

positional arguments:
  vcf                   Path to the (iPyrad) VCF input file.

optional arguments:
  -o OUT_FILE, --out_file OUT_FILE
                        Distances output Phylip file. By default, the distances will be written
                        to a file with the same name as 'vcf' and file extension '.dist'.
  -n N_COL_SKIP, --n_col_skip N_COL_SKIP
                        Number of columns to skip before sample-specific variant data is encountered
                        in VCF file. For iPyrad's VCF output this is 9.
  -e ENTRY_IDX, --entry_idx ENTRY_IDX
                        Index of the variant counts within the variant columns, i.e., number of ':'
                        before variant count information is encountered. For iPyrad's VCF output this
                        is 2. An "entry" with entry_idx=2 might look like this '1/0:295:0,143,0,152'.
  -v, --version         show program's version number and exit
  -h, --help            Show this help message and exit.
```
