#Assume that chemical compund structures are stored in 'molecules.sdf'
from Chem_Faiss import pipeline
import Chem_Faiss
searcher = pipeline('molecules.sdf')

# s is the SMILES string we are looking for in the index; I is the array of matched indexes.
s = 'CCC(CC)COC(=O)[C@H](C)N[P@](=O)(OC[C@H]1O[C@@](C#N)(c2ccc3c(N)ncnn23)[C@H](O)[C@@H]1O)Oc1ccccc1'
I = searcher.make_query_smiles(q = s, k = 2)

# This code block loads the sdf file into mols, and prints the smiles string corresponding to the indexes in I
mols = Chem_Faiss.load_sdf('molecules.sdf')
smiles = Chem_Faiss.idx_to_smiles(mols, I)
print(smiles)

# This code block performs a bulk query from a list of SMILES strings and prints the SMILES string corresponding to the matched indexes in I.
smiles = ['CCC(CC)COC(=O)[C@H](C)N[P@](=O)(OC[C@H]1O[C@@](C#N)(c2ccc3c(N)ncnn23)[C@H](O)[C@@H]1O)Oc1ccccc1',
  'Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)OP(=O)(O)OC[C@H]2O[C@@H](n3cnc4c(N)ncnc43)C(O)C2O)C(O)C1O',
  'NC(=O)c1ccc[n+]([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](O)[C@@H]3O)[C@@H](O)[C@H]2O)c1',
  'NC(=O)C1=CN([C@@H]2O[C@H](CO[P@](=O)(O)O[P@@](=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](O)[C@@H]3O)[C@@H](O)[C@H]2O)C=CC1']
I = searcher.bulk_query_smiles(smiles_list=smiles, k = 1)
smiles = Chem_Faiss.idx_to_smiles(mols, I)

# This code block saves the old pipeline to the disk in the directory: 'sample', initialize a new search pipeline and loads pipeline back from 'sample'.
searcher.save_pipeline('sample')
new_searcher = pipeline()
new_searcher.load_pipeline('sample')
