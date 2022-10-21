import argparse
import fileinput

from gc_meox_src import remove_derivatization_groups, add_derivatization_groups
from concurrent.futures import ProcessPoolExecutor
from rdkit import Chem


def process_one_mol(n_mol):
    mol, n = n_mol
    return (
        mol[0],
        Chem.MolToSmiles(remove_derivatization_groups(mol[1])),
        {Chem.MolToSmiles(add_derivatization_groups(mol[1])) for _ in range(n)}
    )


def parse_args(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument('-n', '--ncpu', type=int, action='store', help='# of cores to use', default=1)
    parser.add_argument('-r', '--repeat', type=int, action='store',
                        help='# of repeated attempts to derivatize (may return different results)', default=42)
    parser.add_argument('-k', '--keep', action='store_true',
                        help='keep input and stripped derivatization SMILES in output', default=True)
    parser.add_argument('-f', '--flat', type=str, action='store', help='flat output file, one SMILES per line')
    parser.add_argument('-t', '--tsv', type=str, action='store',
                        help='structured output tsv file (original, stripped derivatization, added derivatizations')
    parser.add_argument('infiles', nargs='+', type=str, action='store', help='input files')

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)

    insmi = list(filter(lambda p: p[1],
                        [(line.rstrip(), Chem.MolFromSmiles(line)) for line in fileinput.input(files=args.infiles)]))
    n_mols = list(zip(insmi, [args.repeat] * len(insmi)))

    with ProcessPoolExecutor(max_workers=args.ncpu) as executor:
        out = executor.map(process_one_mol, n_mols)

    if args.flat:
        with open(args.flat, "w") as flat:
            if args.keep:
                for orig, removed, added in out:
                    for one in {orig, removed, *added}:
                        flat.write(one + "\n")
            else:
                for orig, removed, added in out:
                    flat.write("\n".join(added) + "\n")

    if args.tsv:
        with open(args.tsv, "w") as tsv:
            tsv.write("orig\tderiv. removed\tderiv. added ...\n")
            for orig, removed, added in out:
                tsv.write("\t".join([orig, removed, *added]) + "\n")


if __name__ == '__main__':
    main()
