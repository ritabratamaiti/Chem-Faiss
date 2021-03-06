from flask import Flask, jsonify, request
from Chem_Faiss import pipeline
import Chem_Faiss


app = Flask(__name__)
searcher = pipeline()
searcher.load_pipeline('sample')
mols = Chem_Faiss.load_sdf('molecules.sdf')


@app.route('/query', methods = ['GET','POST'])
def query():
    d = request.get_json()
    s = d['SMILES']
    print(s)
    print(request.form)
    I = searcher.make_query_smiles(q = s, k = 2)
    smiles = Chem_Faiss.idx_to_smiles(mols, I)
    return jsonify({
    'result' : smiles
    })

app.run()
