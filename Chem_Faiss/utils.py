from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Draw
from sklearn.cluster import KMeans
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.AllChem import *
from tqdm import tqdm
import numpy as np
import faiss
import json
import os
from .utils import *
RDLogger.DisableLog('rdApp.*')

def load_sdf(sdf):
  mols = Chem.SDMolSupplier(sdf)
  n = len(mols)
  mols = [mol for mol in mols if mol is not None]
  print(str(len(mols))+' loaded, '+str(n-len(mols))+' entries have been discarded, as they could not be read by RDKit.')
  return mols

def load_smiles(smiles):
  mols = []
  for s in smiles:
    mols.append(Chem.MolFromSmiles(s))
  n = len(mols)
  mols = [mol for mol in mols if mol is not None]
  print(str(len(mols))+' loaded.\nSince '+str(n-len(mols))+' entries could not be read by RDKit they have been discarded.')
  return mols

def make_fps(mols, fingerprinter = AllChem.GetMACCSKeysFingerprint):
  fps = []
  print('Fingerprinting.....')
  for mol in mols:
      bv = fingerprinter(mol)
      fp = np.zeros(len(bv))
      DataStructs.ConvertToNumpyArray(bv, fp)
      fps.append(fp)
  fps = np.asarray(fps, dtype=np.float32)
  print('Done!')
  return fps

def make_index(fps):
  index = faiss.IndexFlatL2(fps.shape[1])
  index.add(fps)
  print(str(index.ntotal)+' fingerprints added to index.')
  return index

def search(index, mol, k, fingerprinter = AllChem.GetMACCSKeysFingerprint):
  bv = fingerprinter(mol)
  fp = np.zeros(len(bv))
  DataStructs.ConvertToNumpyArray(bv, fp)
  fp = np.asarray([fp], dtype=np.float32)
  D, I = index.search(fp, k)
  return D,I

def bulk_search(index, mols, k, fingerprinter = AllChem.GetMACCSKeysFingerprint):
  fps = make_fps(mols, fingerprinter)
  D, I = index.search(fps, k)
  return D,I

def idx_to_smiles(mols, I):
  temp = []
  for i in range(I.shape[0]):
    temp.append([Chem.MolToSmiles(mols[idx]) for idx in I[i]])
  return temp

def idx_to_mols(mols, I):
  temp = []
  for i in range(I.shape[0]):
    temp.append([mols[idx] for idx in I[i]])
  return temp

def save_index(index, path):
  faiss.write_index(index, path)
  print('Saved to '+path)
def load_index(path):
  index = faiss.read_index(path)
  return index
