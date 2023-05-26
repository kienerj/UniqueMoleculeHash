# Unique Molecule Hash

The intent of `unique_molecule_hash` is to create a unique hash for rdkit molecule instances that can be used to compare them if they match. This should ultimately include every valid rdkit molecule and all it's chemically relevant features including queries.

Existing hashing mechanisms do not always fully suit the needs. inchi does not capture enhanced stereo or dative bonds while RDKit registration hash for example does not capture query features and has other issues that make it not usable for my use-case namely tautomer insensitive hash losing double bond stereo on molecules without tautomerism. 

Unique Molecule Hashes goal is to solve this problem so one can easily compare molecules without having to think which is the right hash to choose. It should "just work". of course this is highly opinionated what should be part of the hash and what not.

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

- **Tautomer-insensitivity** depends on RDKits `TautomerEnumerator` which has several known issues. This can lead to a different hash for the same molecule but with a different input (different kekulization in input format). See [issue 5937](https://github.com/rdkit/rdkit/issues/5937).
- Much slower to generate than InChI
- Same query (SMARTS) can lead to different hash probably due to SMILES canonicalization treating [C] and [C&R1] the same so the order is not guaranteed and depends on input order: `[C]!@;:C-,=C(-[C&R1])-C` vs `"C-C(-[C&R1])-,=C!@;:[C]`. However `[C]!@;:C-,=C(-[C&R1])-C` and `[C]:;!@C-,=C(-[R1&C])-C` do return the same hash so the order of the query does not matter

## Features / Options

The normal way of using `unique_molecule_hash` is just calling `get_standard_hash(mol)` and it will generated a standard hash similar to "standard inChI".

For further configurability there is also the `get_hash` method with several options:


- enumerator: tautomer enumerator to use

    The default one uses default settings.

- tautomer_sensitive: False by default. Whether hash is sensitive to different tautomers of the same molecule

- cx_smiles_fields: flags for cxsmiles creation

    This of course as a big impact on the hash. Coordinates for example are of course excluded in the standard hash.

- normalize_dative_bonds: converts potential dative bonds drawn as single bonds to dative bonds

    Chemists sometimes draw dative bonds as single bonds (leading to invalid valence) or don't draw them at all. The normalization will convert single bonds into dative bonds if applicable and then as a last step remove all dative bonds. This is to get an equal result for all scenarios: drawing it correctly, drawing single bonds instead of dative bonds or not drawing any bonds at all even if dative bond would be relevant

- include_query_features: if query features should be part of the hash or not

    Regardless if this is set, some optimizations  to queries will be applied in all cases. For example `[R1&C]` in SMARTS gets an atomic number 0 (or * in SMILES) while `[C&R1]` is correctly set to 6 and a C in SMILES. This will get fixed in all cases. But if this option is set to `false`, the `R1` part will be "lost" and not part of the hash.

- hash_size: 64 or 128 bits. By default it's 128 bits with very high uniqueness guarantees.

- seed for xxhash

    This just changes the final hash but not the underlying value used for hashing. Can be used when you want to avoid easy look-up from 3rd parties, kind of an anonymization feature

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

Different enhanced stereochemistry leads to a different standard hash

```python
es1 = Chem.MolFromSmiles("C[C@@H](O)CC |o1:1|")
es2 = Chem.MolFromSmiles("C[C@@H](O)CC |a:1|")
es3 = Chem.MolFromSmiles("C[C@@H](O)CC |&1:1|")
unique_molecule_hash.get_standard_hash(es1)
# 07c147d18dc97fc20ebbd710c40743aa
unique_molecule_hash.get_standard_hash(es2)
# d2ed07bf023d09dfafc8f2c8d866cce6
unique_molecule_hash.get_standard_hash(es3)
# b01debf8c10bac48c1af32e49643210d
```

Note that `C[C@@H](O)CC |&1:1`(both stereoisomers) does not get the same hash as drawing both isomers individually: `C[C@@H](O)CC.C[C@H](O)CC`. While one could argue they have the same chemical intend, for now this is out-of-scope but still a possible interesting future enhancement. 

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

