import fileinput
from os import PathLike
from typing import List, Tuple

from rdkit.Chem import Mol, MolFromSmiles


def read_input_txt(infiles: PathLike) -> List[Tuple[str, Mol]]:
    """
    Read input from txt files with SMILES.

    :param infiles: Path to input file(s) with SMILES. One SMILES per line.

    :return: List of tuples (molecule string from the input file, RDKit molecule object of that molecule)
    """
    return [(line.rstrip(), MolFromSmiles(line)) for line in fileinput.input(files=infiles)]


def write_tab_separated(tsv_path: PathLike, data) -> None:
    """
    Write output to a tab-separated file.

    :param tsv_path: Path to output file.
    :param data: Tuple of (original SMILES, underivatized SMILES, set of derivatized SMILES)
    """
    with open(tsv_path, "w") as tsv:
        tsv.write("orig\tderiv. removed\tderiv. added ...\n")
        for orig, removed, added in data:
            tsv.write("\t".join([orig, removed, *added]) + "\n")


def write_flat(txt_path: PathLike, data, keep: bool = False) -> None:
    """
    Write output to a txt file with one SMILES per line.

    :param txt_path: Path to output file.
    :param data: Tuple of (original SMILES, underivatized SMILES, set of derivatized SMILES)
    :param keep: Whether to write the original and underivatized SMILES to the output.
    """
    with open(txt_path, "w") as flat:
        if keep:
            for orig, removed, added in data:
                for one in {orig, removed, *added}:
                    flat.write(one + "\n")
        else:
            for orig, removed, added in data:
                flat.write("\n".join(added) + "\n")
