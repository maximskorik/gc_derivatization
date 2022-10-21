import fileinput

from os import PathLike
from rdkit.Chem import Mol, MolFromSmiles


def read_input_txt(infiles: PathLike) -> list[Mol]:
    """Read input from txt files with SMILES."""
    return [MolFromSmiles(line.rstrip()) for line in fileinput.input(files=infiles)]
