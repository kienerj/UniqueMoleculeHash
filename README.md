# Unique Molecule Hash

The intent of `unique_molecule_hash` is to create a unique hash for rdkit molecule instances that can be used to compare them if they match. This should ultimately include every valid rdkit molecule and all it's chemically relevant features including queries.

Existing hashing mechanisms do not always fully suit the needs. inchi does not capture enhanced stereo or dative bonds while RDKit registration hash for example does not capture query features and has other issues that make it not usable for my use-case namely tautomer insensitive hash losing double bond stereo on molecules without tautomerism. 

Unique Molecule Hashes goal is to solve this problem so one can easily compare molecules without having to think which is the right hash to choose. It should "just work".

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
  - rdkit>=2022.09.5 
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

- Tautomer-insensitivity depends on RDKits `TautomerEnumerator` which has several known issues. This can lead to a different hash for the same molecule but with a different input (different kekulization in input format). See [issue 5937](https://github.com/rdkit/rdkit/issues/5937).
- Much slower to generate than inchi
- Same query (SMARTS) can lead to different hash probably due to SMILES canonicalization treating [C] and [C&R1] the same so the order is not guaranteed and depends on input order: `[C]!@;:C-,=C(-[C&R1])-C` vs `"C-C(-[C&R1])-,=C!@;:[C]`. However `[C]!@;:C-,=C(-[C&R1])-C` and `[C]:;!@C-,=C(-[R1&C])-C` do return the same hash so the order of the query does not matter

## Features / Options

The normal way of using `unique_molecule_hash` is just calling `get_standard_hash(mol)` and it will generated a standard hash similar to "standard inchi".

For further configurability there is also the `get_hash` method with several options:


- enumerator: tautomer enumerator to use

    The default one uses default settings

- cx_smiles_fields: flags for cxsmiles creation

    This of course as a big impact on the hash. Coordinates for example are of course excluded in the standard hash.

- normalize_dative_bonds: converts potential dative bonds drawn as single bonds to dative bonds

    Chemists sometimes draw dative bonds as single bonds (leading to invalid valence) or don't draw them at all. The normalization will convert single bonds into dative bonds if applicable and then as a last step remove all dative bonds. This is to get an equal result for all scenarios: drawing it correctly, drawing single bonds instead of dative bonds or not drawing any bonds at all even if dative bond would be relevant

- include_query_features: if query features should be part of the hash or not

    Regardless if this is set, some optimizations  to queries will be applied in all cases. For example `[R1&C]` in SMARTS gets an atomic number 0 (or * in SMILES) while `[C&R1]` is correctly set to 6 and a C in SMILES. This will get fixed in all cases. But if this option is set to `false`, the `R1` part will be "lost" and not part of the hash.

- seed for xxhash

    This just changes the final hash but not the underlying value used for hashing. Can be used when you want to avoid easy look-up from 3rd parties, kind of an anonymization feature

## Examples

#### Multiple Components

```python
mx1 = Chem.MolFromSmiles("CCC(O)C.CCN")
mx2 = Chem.MolFromSmiles("NCC.CC(O)CC")
unique_molecule_hash.get_standard_hash(mx1)
# 3b6931ee7a71d858
unique_molecule_hash.get_standard_hash(mx2)
# 3b6931ee7a71d858
```

#### Ordering of queries

```python
m1 = Chem.MolFromSmarts("[C]!@;:C-,=C(-[C&R1])-C")
m2 = Chem.MolFromSmarts("[C]:;!@C-,=C(-[R1&C])-C")
unique_molecule_hash.get_standard_hash(m1)
# 98c149472b0154e2
unique_molecule_hash.get_standard_hash(m2)
# 98c149472b0154e2
```

