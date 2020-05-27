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

'''
Low level utility functions.
NOTE: Some minor refactoring needs to be done to allow custom fingerprinting functions in the pipeline module.
'''

def load_sdf(sdf):
  '''
  Reads an SDF file and converts it to a list of RDKit mols, and returns the list.
  '''
  mols = Chem.SDMolSupplier(sdf)
  n = len(mols)
  mols = [mol for mol in mols if mol is not None]
  print(str(len(mols))+' loaded, '+str(n-len(mols))+' entries have been discarded, as they could not be read by RDKit.')
  return mols

def load_smiles(smiles):
  '''
  Accept a list of SMILES strings and converts it to a list of RDKit mols, and returns the list.
  '''
  mols = []
  for s in smiles:
    mols.append(Chem.MolFromSmiles(s))
  n = len(mols)
  mols = [mol for mol in mols if mol is not None]
  print(str(len(mols))+' loaded.\nSince '+str(n-len(mols))+' entries could not be read by RDKit they have been discarded.')
  return mols

def make_fps(mols, fingerprinter = AllChem.GetMACCSKeysFingerprint):
  '''
  Accepts a list of RDKit Mols, and convert to a list of fingerprints. Can use fingerprinter from rdkit.Chem.AllChem
  '''
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
  '''
  Takes an array of fingerprints and converts it to a FAISS index. 
  '''
  index = faiss.IndexFlatL2(fps.shape[1])
  index.add(fps)
  print(str(index.ntotal)+' fingerprints added to index.')
  return index

def search(index, mol, k, fingerprinter = AllChem.GetMACCSKeysFingerprint):
  '''
  Queries a single RDKit Mol against a faiss index. The same fingerprinting scheme used for building the 
  origninal index must be used.
  '''
  bv = fingerprinter(mol)
  fp = np.zeros(len(bv))
  DataStructs.ConvertToNumpyArray(bv, fp)
  fp = np.asarray([fp], dtype=np.float32)
  D, I = index.search(fp, k)
  return D,I

def bulk_search(index, mols, k, fingerprinter = AllChem.GetMACCSKeysFingerprint):
  '''
  Queries a list of RDKit Mols against a faiss index. The same fingerprinting scheme used for building the 
  origninal index must be used.
  '''
  fps = make_fps(mols, fingerprinter)
  D, I = index.search(fps, k)
  return D,I

def idx_to_smiles(mols, I):
  '''
  Given an array of indexes, and a list of RDKit mols, this function will generate the smiles strings for the RDKit
  Mols at the corresponding index location.
  '''
  temp = []
  for i in range(I.shape[0]):
    temp.append([Chem.MolToSmiles(mols[idx]) for idx in I[i]])
  return temp

def idx_to_mols(mols, I):
  '''
  Given an array of indexes, and a list of RDKit mols, this function will generate the smiles strings for the RDKit
  Mols at the corresponding index location.
  '''
  temp = []
  for i in range(I.shape[0]):
    temp.append([mols[idx] for idx in I[i]])
  return temp

def save_index(index, path):
  '''
  Saves a Faiss index to a path.
  '''
  faiss.write_index(index, path)
  print('Saved to '+path)

def load_index(path):
  '''
  Loads a Faiss index from the path.
  '''
  index = faiss.read_index(path)
  return index
