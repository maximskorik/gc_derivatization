# In silico derivatization

Package to perform in-silico MeOX + TMS derivatization (as described e.g. in https://doi.org/10.1021/acs.analchem.7b01010):

* Metoxymation: ketone R(C=O)R' and aldehyde karboxyl groups are substituted with C=NO[CH3]
* Trimethylsilylation: in -OH, -SH, -NH2, -NHR, =NH, the hydrogen is substituted with -SiMe3

The substitution needn't happen always, their probability currently hardcoded in the package.
Typically, multiple substitution attempts are run on each input molecule, and all distinct results are gathered.

Known limitation is metoxymation on cycles which should be broken. This is not implemented yet.

Package provides functions:
* `is_derivatized()` checks whether the molecule contains MeOX or TMS groups that are likely to be result of derivatization
* `remove_derivatization_groups()` removes the suspected groups, reconstructing the original molecule
* `add_derivatization_groups()` does the substitution above

All the functions can accept either `mol: rdkit.Chem.rdchem.Mol` or `smiles: str` argument. All return `rdkit.Chem.rdchem.Mol`.

The typical useage is wrapped in the `gc-meox-tms.py` driver script.

See also https://github.com/ljocha/gc-derivatization for example use in Jupyter notebook.
