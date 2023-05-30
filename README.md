# Unique Molecule Hash

The intent of `unique_molecule_hash` is to create a unique hash for rdkit molecule instances that can be used to compare them if they match. This should ultimately include every valid rdkit molecule and all it's chemically relevant features including queries.

Existing hashing mechanisms do not always fully suit these requirements. InChI does not capture enhanced stereo or dative bonds while RDKit RegistrationHash for example does not capture query features. On top of that the tautomer insensitive hash [may lose double bond stereo on molecules without tautomerism](https://github.com/rdkit/rdkit/discussions/6318). 

Unique Molecule Hashes goal is to solve this problem so one can compare molecules without having to think which is the right hash to choose. It should "just work". Of course this is highly opinionated what should be part of the hash and what not. Therefore the hash can also be configured with multiple options, making it "non-standard". 

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

## Features / Options

The normal way of using `unique_molecule_hash` is just calling `get_standard_hash(mol)` and it will generated a standard hash similar to "standard inChI". 

For further configurability there are either convenience methods `get_molecule_hash` and `get_quick_hash` for common settings or the `get_hash` method with all available options:


- **enumerator**: tautomer enumerator to use

    The default one uses default settings.

- **tautomer_sensitive**: False by default. Whether hash is sensitive to different tautomers of the same molecule

- **cx_smiles_fields**: flags for cxsmiles creation

    This of course as a big impact on the hash. Coordinates for example are of course excluded in the standard hash.

- **normalize_dative_bonds**: converts potential dative bonds drawn as single bonds to dative bonds

    Chemists sometimes draw dative bonds as single bonds (leading to invalid valence) or don't draw them at all. The normalization will convert single bonds into dative bonds if applicable and then as a last step remove all dative bonds. This is to get an equal result for all scenarios: drawing it correctly, drawing single bonds instead of dative bonds or not drawing any bonds at all even if dative bond would be relevant

- **include_query_features**: if query features should be part of the hash or not

    This will also disable optimizations  to queries. For example `[R1&C]` in SMARTS gets an atomic number 0 (or * in SMILES) while `[C&R1]` is correctly set to 6 and a C in SMILES. If `include_query_features` is False, these optimizations will be omitted and the same query in a different order will get a different hash.

- **hash_size**: 64 or 128 bits. By default it's 128 bits with very high uniqueness guarantees.

- **seed** for xxhash

    This just changes the final hash but not the underlying value used for hashing. Can be used when you want to avoid easy look-up from 3rd parties, kind of an anonymization feature

`get_molecule_hash` will by default disable query features and dative bond normalization. This would be a convenience method for general small molecules. The advantage of disabling the features is less processing and hence slightly faster hash generation.

`get_quick_hash` disables all "advanced options" for performance. It should only be applied to pre-cleaned / pre-standardized data sets. It is tautomer sensitive, ignores enhanced stereo and query features and does not normalize dative bonds.

## Caveats

- Tautomer-insensitivity depends on RDKits `TautomerEnumerator` which has several known issues. This can lead to a different hash for the same molecule but with a different input. See [issue 5937](https://github.com/rdkit/rdkit/issues/5937). The issue has been addressed with a workaround but it is suggested you test feasibility on your set of molecules
- Much slower to generate than InChI (also due to workarounds like above)
- The has internally relies on SMILES canonicalization which ignores additional context of atoms like query features and therefore the atom order is not guaranteed and can be different for a different input order.  [See RDKit Issue #6401](https://github.com/rdkit/rdkit/issues/6401). This affects SMARTS with a different input order or multi-center attachments (see the above link for examples.). It does not affect normal small molecules.

## Examples

#### Multiple Components

Multiple components in different order and with different SMILES input

```python
mx1 = Chem.MolFromSmiles("CCC(O)C.CCN")
mx2 = Chem.MolFromSmiles("NCC.CC(O)CC")
unique_molecule_hash.get_standard_hash(mx1)
# 7539c7e1a223bf662968da1579538a03
unique_molecule_hash.get_standard_hash(mx2)
# 7539c7e1a223bf662968da1579538a03
```

#### Ordering of queries

Different order but same meaning of query, bond or atom, leads to same hash

```python
m1 = Chem.MolFromSmarts("[C]!@;:C-,=C(-[C&R1])-C")
m2 = Chem.MolFromSmarts("[C]:;!@C-,=C(-[R1&C])-C")
unique_molecule_hash.get_standard_hash(m1)
# 044da709a16297d9563034b7fffbf134
unique_molecule_hash.get_standard_hash(m2)
# 044da709a16297d9563034b7fffbf134
```

#### Dative Bonds

Molecule drawn with dative bonds, without dative bonds or with single bonds gives same hash

```python
# Correct with Dative Bonds
db1 = Chem.MolFromSmiles("C1=N2->[Mn]3(<-N(=CCCO3)CC2)OCC1")
# No Dative Bonds
db2 = Chem.MolFromSmiles("C1=NCCN=CCCO[Mn]OCC1")
# With single bonds 
# Needs special treatment due to invalid valence
db3 = Chem.MolFromSmiles("[N]12=CCCO[Mn]31OCCC=[N]3CC2", sanitize=False)
db3.UpdatePropertyCache(strict=False)
Chem.SanitizeMol(db3,Chem.SanitizeFlags.SANITIZE_FINDRADICALS|Chem.SanitizeFlags.SANITIZE_KEKULIZE|Chem.SanitizeFlags.SANITIZE_SETAROMATICITY|Chem.SanitizeFlags.SANITIZE_SETCONJUGATION|Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION|Chem.SanitizeFlags.SANITIZE_SYMMRINGS,catchErrors=True)

unique_molecule_hash.get_standard_hash(db1)
# 071ae8601f82ade92537db23c0ff0027
unique_molecule_hash.get_standard_hash(db2)
# 071ae8601f82ade92537db23c0ff0027
unique_molecule_hash.get_standard_hash(db3)
# 071ae8601f82ade92537db23c0ff0027
```
![Alt text](images/dative_bonds.png?raw=true "Polymers")

#### Enhanced Stereochemistry

Different enhanced stereochemistry (and vs or) leads to a different standard hash. 

Different bond direction but same molecule and same enhanced stereo lead to the same hash:

```python
or1 = Chem.MolFromSmiles("C[C@@H](O)CC |o1:1|")
or2 = Chem.MolFromSmiles("C[C@H](O)CC |o1:1|")
unique_molecule_hash.get_standard_hash(or1)
# 07c147d18dc97fc20ebbd710c40743aa
unique_molecule_hash.get_standard_hash(or2)
# 07c147d18dc97fc20ebbd710c40743aa

and1 = Chem.MolFromSmiles("C[C@@H](O)CC |&1:1|")
and2 = Chem.MolFromSmiles("C[C@H](O)CC |&1:1|")
unique_molecule_hash.get_standard_hash(and1)
# b01debf8c10bac48c1af32e49643210d
unique_molecule_hash.get_standard_hash(and2)
# b01debf8c10bac48c1af32e49643210d

abs1 = Chem.MolFromSmiles("C[C@@H](O)CC |a:1|")
abs2 = Chem.MolFromSmiles("C[C@H](O)CC |a:1|")
unique_molecule_hash.get_standard_hash(abs1)
# d2ed07bf023d09dfafc8f2c8d866cce6
unique_molecule_hash.get_standard_hash(abs2)
# 43a2fcd3d3493dacc43128b5eca0dffd
```

#### Repeating Groups / Polymers

```
pm1 = Chem.MolFromMolBlock(poly1)
pm2 = Chem.MolFromMolBlock(poly2)
Chem.Draw.MolsToGridImage([pm1,pm2])
```

![Alt text](images/repeating_groups.png?raw=true "Polymers")

```
unique_molecule_hash.get_standard_hash(pm1)
# 644dd84f9f500710a63a65974b9b1465
unique_molecule_hash.get_standard_hash(pm2)
# 7f131162eb8fce5c79fe4f3d1e794cf5
```

