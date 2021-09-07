import argparse
import pathlib
from sys import path

import nei_vcf

def build_parser():
    parser = argparse.ArgumentParser(
        prog='''nei_vcf''',
        description='''
            Calculate variant-based Nei distances using a VCF input file.
        ''',
        add_help=False,
    )
    parser.add_argument(
        'vcf',
        help='Path to the (iPyrad) VCF input file.',
        type=str,
    )
    parser.add_argument(
        '-o',
        '--out_file',
        help='''
            Distances output Phylip file. By default, the distances
            will be written to a file with the same name as 'vcf' and
            file extension '.dist'.
        ''',
        type=str,
        default=None,
    )
    parser.add_argument(
        '-n',
        '--n_col_skip',
        help='''
            Number of columns to skip before sample-specific variant
            data is encountered in VCF file. For iPyrad's VCF output
            this is 9.
        ''',
        type=int,
        default=9,
    )
    parser.add_argument(
        '-e',
        '--entry_idx',
        help='''
            Index of the variant counts within the variant columns, i.e.,
            number of ':' before variant count information is encountered.
            For iPyrad's VCF output this is 2. An "entry" with entry_idx=2 
            might look like this '1/0:295:0,143,0,152'.
        ''',
        type=int,
        default=2,
    )
    parser.add_argument(
        '-v',
        '--version',
        action='version',
        version='%(prog)s {version}'.format(version=nei_vcf.__version__)
    )
    parser.add_argument(
        '-h', '--help',
        action='help',
        help='Show this help message and exit.',
    )

    return parser

def main(args=None):
    parser = build_parser()
    args = parser.parse_args(args=args)

    vcf_file = args.vcf
    n_col_skip = args.n_col_skip
    entry_idx = args.entry_idx

    dist_outfile = args.out_file
    if dist_outfile is None:
        dist_outfile = str(pathlib.Path(vcf_file).with_suffix('.dist'))

    distances, sample_names = nei_vcf.nei_vcf(
        vcf_file=vcf_file,
        n_col_skip=n_col_skip,
        entry_index=entry_idx,
    )

    max_name_len = max(len(name) for name in sample_names)
    with open(dist_outfile, 'w') as dist_f:
        dist_f.write(str(distances.shape[0]))
        dist_f.write('\n')
        
        for i in range(distances.shape[0]):
            row_tmp = '\t'.join([str(round(x, 6)).ljust(8) for x in distances[i]])
            s = f'{sample_names[i]:{max_name_len}}\t{row_tmp}\n'
            dist_f.write(s)

    print(f'Distance matrix written to: {dist_outfile}')
