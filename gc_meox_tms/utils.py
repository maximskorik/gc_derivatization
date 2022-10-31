import fileinput
from os import PathLike
from typing import List, Tuple

from rdkit.Chem import Mol, MolFromSmiles


def read_input_txt(infiles: PathLike) -> List[Tuple[str, Mol]]:
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
