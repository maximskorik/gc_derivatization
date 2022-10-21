import pytest

from gc_derivatiozation.utils import read_input_txt
from rdkit.Chem import Mol


@pytest.mark.parametrize("path, num_molecules", [
    ("data/acidic_protons.txt", 2),
    ("data/alcohols.txt", 2),
    ("data/ketones.txt", 5)])
def test_reading_input_from_txt(path, num_molecules):
    """Test reading input from txt files."""
    molecules = read_input_txt(path)

    assert len(molecules) == num_molecules
    assert all(isinstance(mol, Mol) for mol in molecules)
