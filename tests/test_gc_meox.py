import pytest

from gc_meox_src import add_derivatization_groups, is_derivatized, remove_derivatization_groups
from rdkit import Chem


FLAKY_RERUNS = 6


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
    """Return a tuple of (smiles, boolean indicating if the molecule is MeOX or TMS derivatized)."""
    smiles, _is_derivatized = request.param
    return smiles, _is_derivatized


@pytest.fixture(params=[
    ("CC(=O)N([Si](C)(C)C)[Si](C)(C)C", "CC(=O)N[Si](C)(C)C", "CC(N)=O"),
    ("C[Si](C)(C)OC1=CC=CC=C1", None, "OC1=CC=CC=C1"),
    ("C[Si](C)(C)OC1=CC=C(O[Si](C)(C)C)C=C1", "C[Si](C)(C)OC1=CC=C(O)C=C1", "OC1=CC=C(O)C=C1"),
    ("CCO[Si](C)(C)C", None, "CCO"),
    ("CC(=O)O[Si](C)(C)C", None, "CC(=O)O"),
    ("CCCS[Si](C)(C)C", None, "CCCS"),
    ("CCC(C)=NOC", None, "CCC(C)=O"),
    ("CC=NOC", None, "CC=O")
])
def derivatization_groups_data(request):
    """Return a tuple of (smiles of a derivatized molecule, smiles of this molecule with different degree of conversion,
    smiles of the original non-derivatized molecule)."""
    derivatized, alternative, original = request.param
    return derivatized, alternative, original


def test_is_derivatized_from_smiles(is_derivatized_data):
    """Test if the is_derivatized function works with SMILES."""
    smiles, expected = is_derivatized_data
    actual = is_derivatized(smiles=smiles)

    assert actual == expected


def test_is_derivatized_from_mol(is_derivatized_data):
    """Test if the is_derivatized function works with RDKit molecules."""
    smiles, expected = is_derivatized_data
    mol = Chem.MolFromSmiles(smiles)
    actual = is_derivatized(mol=mol)

    assert actual == expected


def test_remove_derivatization_groups_from_smiles(derivatization_groups_data):
    """Test if the remove_derivatization_groups function works with SMILES."""
    smiles, _, expected = derivatization_groups_data
    actual = remove_derivatization_groups(smiles=smiles)
    actual_smiles = Chem.MolToSmiles(actual, kekuleSmiles=True)

    assert actual_smiles == expected


def test_remove_derivatization_groups_from_mol(derivatization_groups_data):
    """Test if the remove_derivatization_groups function works with RDKit molecules."""
    smiles, _, expected = derivatization_groups_data
    mol = Chem.MolFromSmiles(smiles)
    actual = remove_derivatization_groups(mol=mol)
    actual_smiles = Chem.MolToSmiles(actual, kekuleSmiles=True)

    assert actual_smiles == expected


@pytest.mark.flaky(reruns=FLAKY_RERUNS)
def test_add_derivatization_groups_from_smiles(derivatization_groups_data):
    """Test if the add_derivatization_groups function works with SMILES. The test will run FLAKY_RERUNS times or until
    success due to non-deterministic nature of add_derivatization_groups."""
    expected, alternative, original = derivatization_groups_data
    derivatized = add_derivatization_groups(smiles=original)
    derivatized_smiles = Chem.MolToSmiles(derivatized, kekuleSmiles=True)

    assert derivatized_smiles in [expected, alternative]


@pytest.mark.flaky(reruns=FLAKY_RERUNS)
def test_add_derivatization_groups_from_mol(derivatization_groups_data):
    """Test if the add_derivatization_groups function works with RDKit molecules. The test will run FLAKY_RERUNS times
    or until success due to non-deterministic nature of add_derivatization_groups."""
    expected, alternative, original = derivatization_groups_data
    mol = Chem.MolFromSmiles(original)
    derivatized = add_derivatization_groups(mol=mol)
    derivatized_smiles = Chem.MolToSmiles(derivatized, kekuleSmiles=True)

    assert derivatized_smiles in [expected, alternative]
