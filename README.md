# Unique Molecule Hash

The intend of `unique_molecule_hash` is to create a unique hash for rdkit molecule instances that can be used to compare them if they match. This should ultimately include every valid rdkit molecule and all it's chemically relevant features including queries.

Existing hashing mechanisms do not always fully suit the needs. inchi does not capture enhanced stereo or dative bonds while RDKit registration hash for example does not capture query features and has other issues that make it not usable for my use-case namely tautomer insensitive hash losing double bond stereo on molecules without tautomerism. 

Unique Molecule Hashes goal is to solve this problem so one can easily compare molecules without having to think which is the right hash to choose. 

This is a work in progress. Use at own risk. Hashes may still change with every new version and there are known limitations.

## Installation

The suggested approach to try it out is to create a new conda environment from an environment.yml:

```yaml
name: UniqueMoleculeHash
channels:  
  - conda-forge 
  - defaults   
dependencies:
  - python>=3.10  
  - rdkit>=2022.09.3 
  - xxhash>=3.2.0
```

```bash
conda env create -f environment.yml
```

Then clone this repository and install it in [development mode](https://packaging.python.org/tutorials/installing-packages/#installing-from-a-local-src-tree) into this new environment using pip:

```bash
python -m pip install -e c:\path\to\UniqueMoleculeHash
```

If you clone the repo with git (vs downloading), this has the advantage that you can use git to pull new commits. 
You can then immediately use the new version without any further changes.

## Caveats

- Tautomer-insensitivity depends on RDKits `TautomerEnumerator` which has several known issues. This can lead to a different hash for the same molecule but with a different input (different kekulization in input format)
- Atom-Query features completely missing for now
- Much slower to generate than inchi
