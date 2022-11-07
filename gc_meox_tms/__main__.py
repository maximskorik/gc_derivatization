import argparse
import sys
from concurrent.futures import ProcessPoolExecutor
from functools import partial

from .derivatization import process_one_mol
from .utils import read_input_txt, write_flat, write_tab_separated


def parse_arguments(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument('-n', '--ncpu', type=int, action='store', help='# of cores to use', default=1)
    parser.add_argument('-r', '--repeat', type=int, action='store',
                        help='# of repeated attempts to derivatize (may return different results)', default=42)
    parser.add_argument('-k', '--keep', action='store_true',
                        help='keep input and stripped derivatization SMILES in output', default=False)
    parser.add_argument('-f', '--flat', type=str, action='store', help='flat output file, one SMILES per line')
    parser.add_argument('-t', '--tsv', type=str, action='store',
                        help='structured output tsv file (original, stripped derivatization, added derivatizations')
    parser.add_argument('infiles', nargs='+', type=str, action='store', help='input files')

    return parser.parse_args(argv)


def main(argv):
    args = parse_arguments(argv)
    input_molecules = read_input_txt(args.infiles)

    process_one_mol_with_repeats = partial(process_one_mol, repeats=args.repeat)
    with ProcessPoolExecutor(max_workers=args.ncpu) as executor:
        data = executor.map(process_one_mol_with_repeats, input_molecules)

    if args.flat:
        write_flat(args.flat, data, args.keep)
    if args.tsv:
        write_tab_separated(args.tsv, data)

    return 0


if __name__ == '__main__':
    main(argv=sys.argv[1:])
