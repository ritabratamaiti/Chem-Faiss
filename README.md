# Chem Faiss

<p align="center">
  <img src="https://github.com/ritabratamaiti/Chem_Faiss/blob/master/Chem_Faiss_logo.png"/>
</p>

This project utilises the vector similarity search functionality from [Faiss](http://https://github.com/facebookresearch/faiss "Faiss"), in conjunction with chemical fingerprinting to build a scalable similarity search architecture for compounds/molecules. This repository also contains an index created from the [ChEMBL SARS-CoV-2 compound data](http://https://www.ebi.ac.uk/chembl/g/#browse/compounds/filter/_metadata.compound_records.src_id%3A52 "ChEMBL SARS-CoV-2 compund data"), available under the `sample/` directory of this repository, which you can load and use via the pipeline object. The [Chem_Faiss_flask.py](https://github.com/ritabratamaiti/Chem_Faiss/blob/master/Chem_Faiss_flask.py) demo uses the index from `sample/`. The data used to create the index is stored as `molecules.sdf`
**Typical use cases include drug design, finding structural matches from a large dataset and drug-repurposement.**
Chem Faiss is an experimental tool developed at [AIDDT.de](http://aiddt.de "AIDDT.de"), and released under the MIT License.

# Quickstart

[Start using Faiss Chem in a Google Colab notebook! Just click the Open in Colab button at the top of the notebook.](https://github.com/ritabratamaiti/Chem-Faiss/blob/master/Chem_Faiss_Colab.ipynb)

# Requirements

*Note: Pick an installation which best suits your needs and runs on your system.*
**Chem Faiss requires:** 
- [Faiss](https://github.com/facebookresearch/faiss "Faiss")
- [RDKit](https://www.rdkit.org/ "RDKit")

Chem Faiss will help you set up the search architecture quickly. I highly recommend going through Faiss as well as RDKit's documentation if you want to extend or modify Chem Faiss to better suit your needs.
You also need [Flask](https://flask.palletsprojects.com/en/1.1.x/ "Flask") if you want to run Chem Faiss as a web service.

# How to use

Chem Faiss allows you to build a FAISS index by simple passing an SDF file or a list of SMILES strings to the Pipeline object. If you have chemical data in a different format, you can use RDKit to read the data as RDKit Mols and then pass this list to the pipeline object.

**Making the index:**

	# Assume that chemical compund structures are stored in 'molecules.sdf'
    from Chem_Faiss import pipeline
    import Chem_Faiss
    searcher = pipeline('molecules.sdf')
**Now you can use the pipeline object (searcher) to query the index like so:**


    # s is the SMILES string we are looking for in the index
    s = 'CCC(CC)COC(=O)[C@H](C)N[P@](=O)(OC[C@H]1O[C@@](C#N)(c2ccc3c(N)ncnn23)[C@H](O)[C@@H]1O)Oc1ccccc1'
    I = searcher.make_query_smiles(q = s, k = 2)
To specify the number of results to be returned, use the argument k. Since k=2 here, this query will return 2 of the most similar molecules from the index. Now that we have the required indexes, we can obtain the actual RDKit Mols/SMILES strings corresponding to this index.


    mols = Chem_Faiss.load_sdf('molecules.sdf')
    smiles = Chem_Faiss.idx_to_smiles(mols, I)
Alternatively, using `Chem_Faiss.idx_to_mols(mols, I)` will generate a list of RDKit mols corresponding to the indexes in `I`.

**You can also save and load the pipeline.** 

    searcher.save_pipeline('sample')
    new_searcher = pipeline()
    new_searcher.load_pipeline('sample')
*Note: Currently, we recommend using the default fingerprinter (rdkit.AllChem.GetMACCSKeysFingerprint) or a fingerprinter provided by RDKit in rdkit.AllChem module, as it leads to greater stability while saving and reloading.*

You can look through [Chem_Faiss_example.py](https://github.com/ritabratamaiti/Chem_Faiss/blob/master/Chem_Faiss_example.py) for a more thorough demo of using Chem Faiss pipeline object. Alternatively, you can use the base utility functions of Chem Faiss to generate your index. These functions will give you greater control over index creation and usage. An example of how to do so can be found in [index_from_scratch.py](https://github.com/ritabratamaiti/Chem_Faiss/blob/master/index_from_scratch.py)
**Note: Some of the functionality in Chem Faiss may not be documented. However, this is an important part of the roadmap.**

## Use Chem Faiss as a web service with Flask!

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

To make a query, send a JSON POST request to `server/query` like so:


    {
    	'SMILES' : 'CCC(CC)COC(=O)[C@H](C)N[P@](=O)(OC[C@H]1O[C@@](C#N)(c2ccc3c(N)ncnn23)[C@H](O)[C@@H]1O)Oc1ccccc1'
    }

## Reference

If you use Chem Faiss in your research, this is the reference to cite:

[![DOI](https://zenodo.org/badge/267277104.svg)](https://zenodo.org/badge/latestdoi/267277104)


## Acknowledgements
- [Faiss](https://github.com/facebookresearch/faiss "Faiss")
- [RDKit](https://www.rdkit.org/ "RDKit")
- [AIDDT](aiddt.de)

