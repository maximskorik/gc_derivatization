import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
from concurrent.futures import ProcessPoolExecutor
from gc_meox_src import is_derivatized, remove_derivatization_groups, add_derivatization_groups
import multiprocessing
import random
import py3Dmol

# define number of cpus to run 
def num_cpu():
    random.seed(42)
    cpus = multiprocessing.cpu_count()
#    print('# cpus (including HT, typically): ', cpus)
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


# reading input file
def read_input_file(filename):
    smi_file = filename
    with open(smi_file) as f:
        mols = list(filter(lambda p: p[1], [ (smi.rstrip(), Chem.MolFromSmiles(smi)) for smi in f ]))
        return mols


# Atcual applcation of the code
def process_one_mol(mol):
    return (
        mol[0],
        Chem.MolToSmiles(remove_derivatization_groups(mol[1])),
        { Chem.MolToSmiles(add_derivatization_groups(mol[1])) for _ in range(42) }
        )


def process_pool_exe(filename):
    with ProcessPoolExecutor(max_workers=num_cpu()) as executor:
        out = executor.map(process_one_mol,read_input_file(filename))
        return out
    

#Writing outputs
def write_derivs_struct(filename):
    with open('derivs_struct.tsv','w') as tsv:
        tsv.write("orig\tderiv. removed\tderiv. added ...\n")
        for orig,removed,added in process_pool_exe(filename):
            tsv.write("\t".join([orig,removed,*added]) + "\n")


def write_derivs_flat(filename):
    with open('derivs_flat.txt','w') as flat:
        for orig,removed,added in process_pool_exe(filename):
            for one in { orig, removed, *added }:
                flat.write(one + "\n")


listarg = argparse.ArgumentParser()
listarg.add_argument('--filename', type=str)
args = listarg.parse_args()
filename = args.filename

if __name__ == "__main__":
    write_derivs_struct(filename) 
    write_derivs_flat(filename)