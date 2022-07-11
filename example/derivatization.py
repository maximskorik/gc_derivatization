from rdkit import Chem
from rdkit.Chem import AllChem
from copy import deepcopy
from concurrent.futures import ProcessPoolExecutor
from gc_meox_src import is_derivatized, remove_derivatization_groups, add_derivatization_groups
import multiprocessing
import random
import py3Dmol

random.seed(42)
cpus = multiprocessing.cpu_count()
print('# cpus (including HT, typically): ', cpus)
cpus //= 2

# 3D rendering
def draw3d(m,dimensions=(500,300),p=None):
    AllChem.EmbedMultipleConfs(m, clearConfs=True, numConfs=50)
    opt = AllChem.MMFFOptimizeMoleculeConfs(m)
    conf = min(range(len(opt)),key = lambda x: opt[x][1] if opt[x][0] == 0 else float("inf") )
    
    mb = Chem.MolToMolBlock(m,confId=conf)
    if p is None:
        p = py3Dmol.view(width=dimensions[0],height=dimensions[1])
    p.removeAllModels()
    p.addModel(mb, 'sdf')
    p.setStyle({'stick': {}})
    p.setBackgroundColor('0xeeeeee')
    p.zoomTo()
    return p.show()


# TESTING imported modules
# is_derivatized
for s in ['CCC(=NOC)C', 'CCC=NOC', 'C=NOC', 'CSi(C)(C)C']:
    print(s,is_derivatized(smiles='CCC(=NOC)C'))

# remove_derivatization 
remove_derivatization_groups(smiles='CCC(=N)C')

m=Chem.MolFromSmiles('CCC=NOC')
remove_derivatization_groups(m)
draw3d(m)
remove_derivatization_groups(smiles='C[Si](C)(C)OCCCO[Si](C)(C)C')

m=remove_derivatization_groups(smiles='CON=CC(O)C=NOC')

# add_derivatization
add_derivatization_groups(m)


# reading input file
smi_file='NIST_195_200.txt'
with open(smi_file) as f:
    mols = list(filter(lambda p: p[1], [ (smi.rstrip(), Chem.MolFromSmiles(smi)) for smi in f ]))


#Essential statistics

SiMe1=Chem.MolFromSmarts('[Si][CH3]')
SiMe2=Chem.MolFromSmarts('[Si]([CH3])[CH3]')
SiMe3=Chem.MolFromSmarts('[Si]([CH3])([CH3])[CH3]')
ONSSi=Chem.MolFromSmarts('[O,N,S][Si]([CH3])([CH3])[CH3]')

print('# total',len(mols))
with_sime1 = list(filter(lambda m: m[1].HasSubstructMatch(SiMe1),mols))
print("# with SiMe:", len(with_sime1))
with_sime2 = list(filter(lambda m: m[1].HasSubstructMatch(SiMe2),mols))
print("# with SiMe2:", len(with_sime2))
with_sime3 = list(filter(lambda m: m[1].HasSubstructMatch(SiMe3),mols))
print("# with SiMe3:", len(with_sime3))
with_onssi = list(filter(lambda m: m[1].HasSubstructMatch(ONSSi),mols))
print("# with ONSSi:", len(with_onssi))

MeOX=Chem.MolFromSmarts('C=NO[CH3]')
with_meox = list(filter(lambda m: m[1].HasSubstructMatch(MeOX),mols))
print("# with MeOX:", len(with_meox))




# Atcual applcation of the code
%%time
def process_one_mol(mol):
    return (
        mol[0],
        Chem.MolToSmiles(remove_derivatization_groups(mol[1])),
        { Chem.MolToSmiles(add_derivatization_groups(mol[1])) for _ in range(42) }
        )
        
with ProcessPoolExecutor(max_workers=cpus) as executor:
    out = executor.map(process_one_mol,mols)
    
out = list(out)

#Writing outputs
with open('derivs_struct.tsv','w') as tsv:
    tsv.write("orig\tderiv. removed\tderiv. added ...\n")
    for orig,removed,added in out:
        tsv.write("\t".join([orig,removed,*added]) + "\n")

with open('derivs_flat.txt','w') as flat:
    for orig,removed,added in out:
        for one in { orig, removed, *added }:
            flat.write(one + "\n")