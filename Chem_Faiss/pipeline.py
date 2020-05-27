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

class pipeline:
  '''
  Creates a pipeline to make an index, search queries, and save indexes to the disk. 
  '''
  def __init__(self, sdf_path = None, smiles_list = None, mols_list = None, fingerprinter = None):
    '''
    Initialises a pipeline and creates an index. Can be configured to use a custom fingerprinter. Accepts SDF files,
    lists of SMILES strings or a list of RDKit Mols to generate the index.
    * NOTE: Typically works with a fingerprinter defined in rdkit.Chem.AllChem
    * NOTE: To be refactored in the future for conciseness.
    '''
    self.index = None
    if fingerprinter is None:
      self.fingerprinter = AllChem.GetMACCSKeysFingerprint
      if sdf_path:
        mols = load_sdf(sdf_path)
        fps = make_fps(mols)
        self.index = make_index(fps)

      elif smiles_list:
        mols = load_smiles(smiles_list)
        fps = make_fps(mols)
        self.index = make_index(fps)

      elif mols_list:
        fps = make_fps(mols)
        self.index = make_index(fps)
    else:
      self.fingerprinter = fingerprinter
      if sdf_path:
        mols = load_sdf(sdf_path)
        fps = make_fps(mols, fingerprinter)
        self.index = make_index(fps)
        return mols
      elif smiles_list:
        mols = load_smiles(smiles_list)
        fps = make_fps(mols, fingerprinter)
        self.index = make_index(fps)
        return mols
      elif mols_list:
        fps = make_fps(mols, fingerprinter)
        self.index = make_index(fps)
  def make_query_smiles(self, q, k = 4):
    '''
    Query the index for a single SMILES string. This will return top k results from the index.
    '''
    mol = Chem.MolFromSmiles(q)
    D, I = search(self.index, mol, k, self.fingerprinter)
    return I
  def make_query_mol(self, mol, k = 4):
    '''
    Query the index for a single RDKit Mol. This will return top k results from the index.
    '''
    D, I = search(self.index, mol, k, self.fingerprinter)
    return I
  def bulk_query_smiles(self, smiles_list, k = 4):
    '''
    Query the index for a list of SMILES string. This will return top k results from the index for each of the list entries.
    '''
    mols = load_smiles(smiles_list)
    D, I = bulk_search(self.index, mols, k, self.fingerprinter)
    return I
  def bulk_query_mols(self, mols, k = 4):
    '''
    Query the index for a list of RDKit Mols. This will return top k results from the index for each of the list entries.
    '''
    D, I = bulk_search(self.index, mols, k, self.fingerprinter)
    return I
  def save_pipeline(self, path):
    '''
    Saves the pipeline in the path directory. The index and a JSON config file will be saved.
    '''
    if not os.path.exists(path):
      os.makedirs(path)
    else:
      print('This directory already exists, please choose a different name!')
    save_index(self.index, path+'/search.index')
    data = {'fingerprinter': self.fingerprinter.__name__}
    with open(path+'/config.json', 'w', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=4)
  def load_pipeline(self, path):
    '''
    Loads the pipeline from the path directory. 
    '''
    self.index = load_index(path+'/search.index')
    with open(path+'/config.json') as f:
        data = json.load(f)
    self.fingerprinter = eval(data['fingerprinter'])
