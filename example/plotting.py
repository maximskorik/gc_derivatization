import py3Dmol

from rdkit.Chem import AllChem, Mol, MolToMolBlock
from typing import Optional


def draw3d(m: Mol, dimensions: tuple[int, int] = (500, 300), p: Optional[py3Dmol.view] = None):
    AllChem.EmbedMultipleConfs(m, clearConfs=True, numConfs=50)
    opt = AllChem.MMFFOptimizeMoleculeConfs(m)
    conf = min(range(len(opt)), key=lambda x: opt[x][1] if opt[x][0] == 0 else float("inf"))
    mb = MolToMolBlock(m, confId=conf)

    if p is None:
        p = py3Dmol.view(width=dimensions[0], height=dimensions[1])

    p.removeAllModels()
    p.addModel(mb, 'sdf')
    p.setStyle({'stick': {}})
    p.setBackgroundColor('0xeeeeee')
    p.zoomTo()
    return p.show()
