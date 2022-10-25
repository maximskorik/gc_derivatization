import pytest

from gc_meox_src import is_derivatized, remove_derivatization_groups
from rdkit import Chem


@pytest.fixture(params=[
    ("CC(=O)N([Si](C)(C)C)[Si](C)(C)C", True),
    ("C[Si](C)(C)OC1=CC=CC=C1", True),
    ("C[Si](C)(C)OC1=CC=C(C=C1)O[Si](C)(C)C", True),
    ("C[Si](C)(C)C1=CC=C(C=C1)[Si](C)(C)C", False),
    ("CCO[Si](C)(C)C", True),
    ("CC(=O)O[Si](C)(C)C", True),
    ("CC(=O)O", False),
    ("CCCS[Si](C)(C)C", True),
    ("CCCS", False),
    ("CCC(=NOC)C", True),
    ("CC=NOC", True),
    ("CCCC(=O)N", False),
    ("CCCC(=O)NCC", False),
    ("CC(=O)NOC", False),
    ("CCC(O)C", False),
    ("CCCC#N", False),
    ("C[N+]#[C-]", False)
])
def is_derivatized_data(request):
    smiles = request.param[0]
    expected = request.param[1]
    return smiles, expected


@pytest.fixture(params=[
    ("CC(=O)N([Si](C)(C)C)[Si](C)(C)C", "CC(N)=O"),
    ("C[Si](C)(C)OC1=CC=CC=C1", "OC1=CC=CC=C1"),
    ("C[Si](C)(C)OC1=CC=C(C=C1)O[Si](C)(C)C", "OC1=CC=C(O)C=C1"),
    ("CCO[Si](C)(C)C", "CCO"),
    ("CC(=O)O[Si](C)(C)C", "CC(=O)O"),
    ("CCCS[Si](C)(C)C", "CCCS"),
    ("CCC(=NOC)C", "CCC(C)=O"),
    ("CC=NOC", "CC=O")
])
def derivatization_groups_data(request):
    smiles = request.param[0]
    expected = request.param[1]
    return smiles, expected


def test_is_derivatized_from_smiles(is_derivatized_data):
    smiles, expected = is_derivatized_data
    actual = is_derivatized(smiles=smiles)

    assert actual == expected


def test_is_derivatized_from_mol(is_derivatized_data):
    smiles, expected = is_derivatized_data
    mol = Chem.MolFromSmiles(smiles)
    actual = is_derivatized(mol=mol)

    assert actual == expected


def test_remove_derivatization_groups_from_smiles(derivatization_groups_data):
    smiles, expected = derivatization_groups_data
    actual = remove_derivatization_groups(smiles=smiles)
    actual_smiles = Chem.MolToSmiles(actual, kekuleSmiles=True)

    assert actual_smiles == expected


def test_remove_derivatization_groups_from_mol(derivatization_groups_data):
    smiles, expected = derivatization_groups_data
    mol = Chem.MolFromSmiles(smiles)
    actual = remove_derivatization_groups(mol=mol)
    actual_smiles = Chem.MolToSmiles(actual, kekuleSmiles=True)

    assert actual_smiles == expected
