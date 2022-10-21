import fileinput

from gc_meox_src import add_derivatization_groups, remove_derivatization_groups
from os import PathLike
from rdkit.Chem import Mol, MolFromSmiles, MolToSmiles


def process_one_mol(n_mol):
    mol, n = n_mol
    return (
        mol[0],
        MolToSmiles(remove_derivatization_groups(mol[1])),
        {MolToSmiles(add_derivatization_groups(mol[1])) for _ in range(n)}
    )


def read_input_txt(infiles: PathLike) -> list[tuple[str, Mol]]:
    """Read input from txt files with SMILES."""
    return [(line.rstrip(), MolFromSmiles(line)) for line in fileinput.input(files=infiles)]


def write_tab_separated(tsv_path: PathLike, data) -> None:
    with open(tsv_path, "w") as tsv:
        tsv.write("orig\tderiv. removed\tderiv. added ...\n")
        for orig, removed, added in data:
            tsv.write("\t".join([orig, removed, *added]) + "\n")


def write_flat(txt_path: PathLike, data, keep: bool = False) -> None:
    with open(txt_path, "w") as flat:
        if keep:
            for orig, removed, added in data:
                for one in {orig, removed, *added}:
                    flat.write(one + "\n")
        else:
            for orig, removed, added in data:
                flat.write("\n".join(added) + "\n")
