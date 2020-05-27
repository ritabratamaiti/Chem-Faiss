from Chem_Faiss import *
# Please look through Chem_Faiss/utils for implementation details.

#This code block prepares the index.
mols = load_sdf('molecules.sdf')
fps = make_fps(mols)
index = make_index(fps)

#This code block searched the index for a query.
s = 'CCC(CC)COC(=O)[C@H](C)N[P@](=O)(OC[C@H]1O[C@@](C#N)(c2ccc3c(N)ncnn23)[C@H](O)[C@@H]1O)Oc1ccccc1'
D, I = search(index, Chem.MolFromSmiles(s), 4)
idx_to_smiles(mols, I)

#Saving and loading an index.
save_index(index, 'sample.index')
new_index = load_index('sample.index')
