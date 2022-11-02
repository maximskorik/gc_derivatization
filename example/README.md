# In silico derivatization example

The example Jupyter notebook `derivatiozation.ipynb` reads a list of SMILES (text file, one molecule per line), and performs the derivatisation, also inspecting its results.

The final outputs are two files:

* `derivs_struct.tsv` with columns (all SMILES):
  * original
  * with derivatization groups stripped
  * column #2 derivatized (multiple times) according to the above rules
* `derivs_flat.txt` -- the above with all the smiles flattened, one per line


