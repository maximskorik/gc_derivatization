import pytest

from gc_meox_tms.utils import read_input_txt
from rdkit.Chem import Mol


@pytest.mark.parametrize("path, smiles", [
    ("data/acidic_protons.txt", ["CC(=O)O", "C(C(C(=O)O)N)S"]),
    ("data/alcohols.txt", ["CCO", "CO"]),
    ("data/ketones.txt", ["CC(=O)C", "CCC(=O)C", "CC(=O)CC(=O)C", "CC1CCCCCCCCCCCCC(=O)C1", "C1CCC(=O)CC1"])
])
def test_reading_input_from_txt(path, smiles):
    """Test reading input from txt files."""
    molecules = read_input_txt(path)
    actual_smiles = [mol[0] for mol in molecules]
    rdkit_molecules = [mol[1] for mol in molecules]

    assert len(molecules) == len(smiles)
    assert actual_smiles == smiles
    assert all(isinstance(mol, Mol) for mol in rdkit_molecules)
