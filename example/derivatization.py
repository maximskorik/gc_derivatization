import argparse
import random
import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem
from concurrent.futures import ProcessPoolExecutor
from gc_meox_src import is_derivatized, remove_derivatization_groups, add_derivatization_groups
import multiprocessing

dict_deriv = {
    "SiMe1" : Chem.MolFromSmarts('[Si][CH3]'),
    "SiMe2" : Chem.MolFromSmarts('[Si]([CH3])[CH3]'),
    "SiMe3" : Chem.MolFromSmarts('[Si]([CH3])([CH3])[CH3]'),
    "ONSSi" : Chem.MolFromSmarts('[O,N,S][Si]([CH3])([CH3])[CH3]'),
    "MeOX"  : Chem.MolFromSmiles('CC1CO1')
}

# define number of cpus to run 
def num_cpu():
    random.seed(42)
    cpus = multiprocessing.cpu_count()
#    print('# cpus (including HT, typically): ', cpus)
    cpus //= 2


# Utility 3D rendering and plot
def set_draw3d(m,dimensions=(500,300),p=None):
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

# works interactively in notebook only   
def plot_3d(key):
    x = list(filter(lambda m: m[1].HasSubstructMatch(dict_deriv[key]),
    read_input_file(filename)))
    set_draw3d(x[0][1])

# reading input file
def read_input_file(filename):
    smi_file = filename
    with open(smi_file) as f:
        mols = list(filter(lambda p: p[1], [ (smi.rstrip(), Chem.MolFromSmiles(smi)) for smi in f ]))
        return mols


# Actual applcation of the code
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
    

# Writing outputs
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


# Essential statistics 
def statistics(filename):
    print('# Total',':',len(read_input_file(filename)))
    for key in dict_deriv:
        x = list(filter(lambda m: m[1].HasSubstructMatch(dict_deriv[key]),
        read_input_file(filename)))
        print('# With',key,':',len(x))


# Check on manual inputs
def manual_check():
    for s in 'CCC(=NOC)C', 'CCC=NOC', 'C=NOC', 'CSi(C)(C)C':
        print(s,is_derivatized(smiles='CCC(=NOC)C'))
        m = remove_derivatization_groups(smiles='CCC(=N)C')
        n = remove_derivatization_groups(smiles='C[Si](C)(C)OCCCO[Si](C)(C)C')
        p = remove_derivatization_groups(smiles='CON=CC(O)C=NOC')
        s = add_derivatization_groups(n)


listarg = argparse.ArgumentParser()
listarg.add_argument('--filename', type=str)
args = listarg.parse_args()
filename = args.filename

if __name__ == "__main__":
    statistics(filename)
    write_derivs_struct(filename) 
    write_derivs_flat(filename)
    plot_3d('SiMe3')
    manual_check()