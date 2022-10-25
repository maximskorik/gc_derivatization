from rdkit import Chem
from rdkit.Chem import AllChem
from copy import deepcopy
import random

tms = '[Si]([CH3])([CH3])[CH3]'

# XXX: ~[O,N,S] would match more than we aim to (-O, -S, -N, =N) but it's unlikely to happen
tms_match = Chem.MolFromSmarts('*~[O,N,S]' + tms)
tms_match0 = Chem.MolFromSmarts('[#0]([CH3])([CH3])[CH3]')

meox_match_co = Chem.MolFromSmarts('C([C,c])([C,c])=NO[CH3]')
meox_match_cho = Chem.MolFromSmarts('[CH]([C,c])=NO[CH3]')
meox_match0 = Chem.MolFromSmarts('[#0]=NO[CH3]')
co = Chem.MolFromSmiles('C=O')


def is_derivatized(mol=None, smiles=None):
    if mol is None:
        mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    return mol.HasSubstructMatch(tms_match) or mol.HasSubstructMatch(meox_match_co) or mol.HasSubstructMatch(
        meox_match_cho)


def remove_derivatization_groups(mol=None, smiles=None):
    if mol is None:
        em = Chem.MolFromSmiles(smiles)
    else:
        em = deepcopy(mol)

    matches = em.GetSubstructMatches(tms_match)
    for ma in matches:
        em.GetAtomWithIdx(ma[2]).SetAtomicNum(0)

    em = AllChem.DeleteSubstructs(em, tms_match0)

    matches = em.GetSubstructMatches(meox_match_co)
    for ma in matches:
        em.GetAtomWithIdx(ma[0]).SetAtomicNum(0)
    matches = em.GetSubstructMatches(meox_match_cho)
    for ma in matches:
        em.GetAtomWithIdx(ma[0]).SetAtomicNum(0)

    em, = AllChem.ReplaceSubstructs(em, meox_match0, co, replaceAll=True)
    Chem.SanitizeMol(em)
    return em


# (match pattern, dummy atom #, probability)
_subs = [
    ('[OH]', [100], [.95]),
    ('[SH]', [101], [.80]),
    # matches also imine
    ('[NH]', [102], [.50]),
    ('[NH2]', [103, 102], [.25, .5]),
    ('C([C,c])([C,c])=O', [104], [.90]),
    ('[CH]=O', [104], [.90]),
]

# (dummy atom #, replacement)
_repls = [
    ('[#100]', 'O' + tms),
    ('[#101]', 'S' + tms),
    ('[#102]', 'N' + tms),
    ('[#103]', f'N({tms}){tms}'),
    ('[#104]=O', 'C=NO[CH3]'),
]

# repls = list(zip(
#    map(lambda n: Chem.MolFromSmarts(f'[#{n}]'),_repls),
#    map(Chem.MolFromSmiles,_repls.values())
# ))

subs = [(Chem.MolFromSmarts(pat), repls, probs) for pat, repls, probs in _subs]
repls = [(Chem.MolFromSmarts(pat), Chem.MolFromSmiles(repl)) for pat, repl in _repls]


def add_derivatization_groups(mol=None, smiles=None):
    if mol is None:
        mol = Chem.MolFromSmiles(smiles)

    em = deepcopy(mol)

    for pat, reps, probs in subs:
        matches = em.GetSubstructMatches(pat)
        #        print(matches)
        for m in matches:
            r = random.random()
            for repl, prob in zip(reps, probs):
                if r < prob:
                    em.GetAtomWithIdx(m[0]).SetAtomicNum(repl)
                    break

    for pat, repl in repls:
        #       print(Chem.MolToSmiles(pat),Chem.MolToSmiles(repl),Chem.MolToSmiles(em))
        em, = AllChem.ReplaceSubstructs(em, pat, repl, replaceAll=True)

    Chem.SanitizeMol(em)
    return em
