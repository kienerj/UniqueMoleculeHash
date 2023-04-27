# UniqueMoleculeHash

`unique_molecule_hash` creates a unique hash for rdkit molecules that can be used to compare them if they match exactly.

## Installation

The suggested approach to try it out is to create a new conda environment from an environment.yml:

```yaml
name: UniqueMoleculeHash
channels:  
  - conda-forge 
  - defaults   
dependencies:
  - python>=3.10  
  - rdkit>=2022.09.1 
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

## Usage Examples


