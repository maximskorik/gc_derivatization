import pytest

from concurrent.futures import ProcessPoolExecutor
from gc_meox_tms.utils import read_input_txt
from pathlib import Path
from rdkit.Chem import Mol, MolFromSmiles
from utils import process_one_mol, write_flat, write_tab_separated


@pytest.fixture(scope="module")
def test_dir(request):
    """Return the directory of the currently running test script."""
    return Path(request.fspath).parent


@pytest.fixture
def data():
    molecules = [(smiles, MolFromSmiles(smiles)) for smiles in [
        "CCC(=NOC)C", "CCC=NOC", "C=NOC", "CC(=O)N([Si](C)(C)C)[Si](C)(C)C"]]
    num_molecules = list(zip(molecules, [1] * len(molecules)))

    with ProcessPoolExecutor(max_workers=2) as executor:
        data = executor.map(process_one_mol, num_molecules)

    yield data


@pytest.mark.parametrize("path, smiles", [
    ("data/acidic_protons.txt", ["CC(=O)O", "C(C(C(=O)O)N)S"]),
    ("data/alcohols.txt", ["CCO", "CO"]),
    ("data/ketones.txt", ["CC(=O)C", "CCC(=O)C", "CC(=O)CC(=O)C", "CC1CCCCCCCCCCCCC(=O)C1", "C1CCC(=O)CC1"])
])
def test_reading_input_from_txt(path, test_dir, smiles):
    """Test reading input from txt files."""
    molecules = read_input_txt(test_dir / path)
    actual_smiles = [mol[0] for mol in molecules]
    rdkit_molecules = [mol[1] for mol in molecules]

    assert len(molecules) == len(smiles)
    assert actual_smiles == smiles
    assert all(isinstance(mol, Mol) for mol in rdkit_molecules)


def test_writing_flat_output(data, tmp_path):
    """Test writing flat output."""
    flat_path = tmp_path / "flat.txt"
    write_flat(flat_path, data, True)

    assert flat_path.exists()


def test_writing_flat_content(data, tmp_path):
    """Test writing flat output content."""
    flat_path = tmp_path / "flat.txt"
    write_flat(flat_path, data, True)

    with open(flat_path, "r") as f:
        lines = f.readlines()

    assert len(lines) == 8


def test_writing_tsv_output(data, tmp_path):
    """Test writing tsv output."""
    tsv_path = tmp_path / "tsv.txt"
    write_tab_separated(tsv_path, data)

    assert tsv_path.exists()


def test_writing_tsv_content(data, tmp_path):
    """Test writing tsv output content."""
    tsv_path = tmp_path / "tsv.txt"
    write_tab_separated(tsv_path, data)

    with open(tsv_path, "r") as f:
        lines = f.readlines()

    assert len(lines) == 5
    assert lines[0] == "orig\tderiv. removed\tderiv. added ...\n"


@pytest.mark.parametrize("smiles, expected", [
    ("CC(N)=O", {"CC(N)=O", "CC(=O)N([Si](C)(C)C)[Si](C)(C)C", "CC(=O)N[Si](C)(C)C"}),
    ("C[Si](C)(C)OC1=CC=C(O)C=C1", {"OC1=CC=C(O)C=C1", "C[Si](C)(C)OC1=CC=C(O[Si](C)(C)C)C=C1", "C[Si](C)(C)OC1=CC=C(O)C=C1"}),
    ("CCC(C)=O", {"CCC(C)=O", "CCC(C)=NOC"}),
    ("CC=NOC", {"CC=O", "CC=NOC"})
])
def test_process_one_mol(smiles, expected):
    """Test processing one molecule."""
    mol = (smiles, MolFromSmiles(smiles))
    n = 40
    actual = process_one_mol((mol, n))
    actual = {actual[0], actual[1], *actual[2]}

    assert actual == expected
