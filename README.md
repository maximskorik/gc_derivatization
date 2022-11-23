# In silico derivatization

## Overview

This package performs in-silico MeOX + TMS derivatization (as described e.g. in https://doi.org/10.1021/acs.analchem.7b01010):

* Methoximation: ketone R(<font color='pink'>C</font>=O)R' and aldehyde (-H<font color='pink'>C</font>=O) karboxyl groups 
are substituted with -<font color='pink'>C</font>=NOCH<sub>3</sub>
* Trimethylsilylation: the acidic hydrogen in -OH, -SH, -COOH, -NH<sub>2</sub>, -NHR, =NH, the hydrogen is substituted with -Si(CH<sub>3</sub>)<sub>3</sub>
The substitution needn't happen always, their probability currently hardcoded in the package.
Typically, multiple substitution attempts are run on each input molecule, and all distinct results are gathered.

Known limitation is methoximation on cycles which should be broken. This is not implemented yet.


## Installation

1. From source by cloning the repository and installing the package with `pip` as follows:
```shell
$ git clone https://github.com/RECETOX/gc-meox-tms.git

# install the package:
$ python -m pip install gc-meox-tms

# if you want to run examples in the Jupyter notebook, install with this command:
$ python -m pip install gc-meox-tms[eda]
```

Other installation ways are not implemented yet.

## Usage

### Command-Line Tool

`gc-meox-tms` can be used as a command line tool to produce all MeOX/TMS derivatives of given compounds. To use it via
the command line you will need one or more `txt` files with chemical compounds represented as SMILES
(one SMILES per line). The tool can output results in flat `txt` format(one compound per line) or tab separated `tsv`
format (all derivatives of a given molecule per line).
```shell
$ python -m gc_meox_tms \
-f <path to write flat txt result> \
-t <path to write tab separated result> \
<paths to input txt files>
```
More parameters can be specified, such as number of cores or repeats. For more information run:
```shell
$ python -m gc_meox_tms --help
```

### Python Package

Package provides functions:
* `is_derivatized()` checks whether the molecule contains MeOX or TMS groups that are likely to be result of derivatization
* `remove_derivatization_groups()` removes the suspected groups, reconstructing the original molecule
* `add_derivatization_groups()` does the substitution above

```python3
from gc_meox_tms import add_derivatization_groups, is_derivatized, remove_derivatization_groups
from rdkit.Chem import MolToSmiles

# Example compounds in SMILES format
compounds = ["CC=O", "CC=NOC", "CCO[Si](C)(C)C"]

# Check derivatization
[is_derivatized(smiles=smiles) for smiles in compounds]
>>> [False, True, True]

# Remove derivatization groups from derivatized molecules
underivatized = [remove_derivatization_groups(smiles=smiles) for smiles in compounds[1:]]
print([MolToSmiles(mol) for mol in underivatized])
>>> ["CC=O", "CCO"]

# Convert molecules back to derivatized forms
rederivatized = [add_derivatization_groups(mol=mol) for mol in underivatized]
print([MolToSmiles(mol) for mol in rederivatized])
>>> ['CC=NOC', 'CCO[Si](C)(C)C']
```
Note that your results may differ from the presented since `add_derivatization_groups` is not deterministic. If you rerun
the function enough times you will get all possible derivatizations. The number of reruns to obtain all possible conformations
is individual for each compound (depends on possible conversion degrees etc.).

See also the Jupyter notebook in `example/` directory for more examples.
