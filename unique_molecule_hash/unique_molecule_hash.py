from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
import xxhash
import rdkit
import re
import logging
import ast


logger = logging.getLogger('unique_molecule_hash')

# Creating this instance is very costly so doing it only once instead of per function call reduces runtime by 5x!
tautomer_enumerator = rdMolStandardize.TautomerEnumerator()
p_atom_type = re.compile(r"AtomType (\d+)")
p_bond_type = re.compile(r"BondOrder (\d+)")


def separate_components(mol: Chem.Mol) -> list:
    """
    Separates the rdkit mol into its individual components and returns them as a list of rdkit RWMols for further
    manipulation.

    This function is required to circumvent a bug in rdkit before 2023.03.1. Chem.GetMolFrags removes enhanced stereo
    from the returned fragments. This reassigns it.

    :param mol: a valid rdkit molecule
    :return: a list of all components/fragments of the input mol
    """
    mapping = []
    frags = []
    logger.debug(f"Separating components for molecule {Chem.MolToSmiles(mol)}")
    components = Chem.GetMolFrags(mol, asMols=True, frags=frags, fragsMolAtomMapping=mapping)

    # tuple to list for later assignments
    components = list(components)
    logger.debug(f"Molecule consists of {len(components)} components.")

    # GetMolFrags prior to 2023.03.1 loses enhanced stereo if there is more than one component (structure) in the
    # molecule. Therefore, the enhanced stereo needs to be recreated and reassigned.
    if rdkit.__version__ < "2023.03.01" and len(components) > 1 and len(mol.GetStereoGroups()) > 0:
        logger.debug(f"Found RDKit version < 2023.03.01. Have to reassign enhanced stereo to components.")
        stereo_groups = [[] for i in range(len(components))]

        for st_grp in mol.GetStereoGroups():
            grp_type = st_grp.GetGroupType()
            atoms = st_grp.GetAtoms()
            # group is  limited to 1 component in case of or / and. but only a single group for abs
            # StereoGroup for absolute stereochemistry is only one group with all atoms from both components
            # Therefore code needs to split them up and create a StereoGroup for each component
            if grp_type == Chem.STEREO_ABSOLUTE:
                atom_ids = [[] for i in range(len(components))]
                for atom in atoms:
                    old_idx = atom.GetIdx()
                    for component_idx, component in enumerate(components):
                        if old_idx in mapping[component_idx]:
                            new_idx = mapping[component_idx].index(old_idx)
                            atom_ids[component_idx].append(new_idx)
                            break
                # Build STEREO_ABSOLUTE group for each component separately
                for component_idx in range(len(components)):
                    if len(atom_ids[component_idx]) > 0:
                        if len(stereo_groups[component_idx]) == 0:
                            rwm = Chem.RWMol(components[component_idx])
                            components[component_idx] = rwm
                        sg = Chem.CreateStereoGroup(st_grp.GetGroupType(), components[component_idx], atom_ids[component_idx])
                        stereo_groups[component_idx].append(sg)
            else:
                a0_idx = atoms[0].GetIdx()
                component_idx = frags[a0_idx]
                atom_ids = []
                for atom in atoms:
                    old_idx = atom.GetIdx()
                    new_idx = mapping[component_idx].index(old_idx)
                    atom_ids.append(new_idx)
                # change component to RWMol and get atoms from this RWMol
                if len(stereo_groups[component_idx]) == 0:
                    rwm = Chem.RWMol(components[component_idx])
                    components[component_idx] = rwm
                sg = Chem.CreateStereoGroup(st_grp.GetGroupType(), components[component_idx], atom_ids)
                stereo_groups[component_idx].append(sg)

        for component_idx, component in enumerate(components):
            if len(stereo_groups[component_idx]) > 0:
                component.SetStereoGroups(stereo_groups[component_idx])
                components[component_idx] = component
                logger.debug(f"Fixed Stereo for component {Chem.MolToSmiles(component)}")
    else:
        components = [Chem.RWMol(comp) for comp in components]
    return components


def get_standard_hash(mol: Chem.Mol):
    """
    Generates a unique hash for this molecule to compare rdkit molecules on tehir same chemical intent.

    :param mol: a valid rdkit molecule
    :return: the stanard unique hash
    """
    return get_hash(mol)


# TODO: include atom queries in hash. unclear how as GetSmarts() is not canonical so the same query with different
# order returns a different SMARTS and the cxsmiles always contains the first atom in the query
# eg [N,O,S] means CXSMiles will contain the "N" while in [O,N,S] it would be the "O".
# Tautomerism: open issues around canonical tautomers not being at all that canonical especially if input is differently
# kekulized (https://github.com/rdkit/rdkit/issues/5937). Use inchi if no query features present?
# idea: make standard hash method and a configurable one to choose tautomer layer (inchi or enumerator), custom layer
# with or without enhanced stereo etc.
def get_hash(mol: Chem.Mol, enumerator=tautomer_enumerator, tautomer_sensitive: bool = False,
             cx_smiles_fields: int = 489, normalize_dative_bonds: bool = True, include_query_features: bool = True,
             seed: int = 0) -> str:
    """
    Creates a hash to compare RDKit molecules. It takes into account enhanced stereo, tautomerism and optionally
    query features.

    By default, an rdMolStandardize.TautomerEnumerator() instance is used. You can pass your own as long
    as it has a "Canonicalize(mol)" method that returns a canonical tautomer (RDKit molecule). This of course impact the
    generated hash.

    By default the hash is tautomer insensitive, meaning a different tautomer of the same molecule will result in the
    same hash. To make the hash tautomer sensitive, set tautomer_sensitive = True.
    Note that the rdMolStandardize.TautomerEnumerator() has bugs in which a different input with different kekulization
    can lead to a different canonical tautomer and hence a different hash!

    normalize_dative_bonds option converts potential dative bonds drawn as single bonds to dative bonds. This is only
    done if invalid valences are found. As a last step all dative bonds will be removed to give the same hash to the
    same "chemical intent" as chemists often don't draw dative bonds at all or draw them as single bonds.

    include_query_features option will include query features as part of the hash. In any case, query atoms will be
    included as "any atom" (*) but with this option the actual query will be part of the hash as well.

    :param mol: a valid rdkit molecule
    :param enumerator: tautomer enumerator to use
    :param tautomer_sensitive: if the hash is sensitive to different tautomers of the same molecule or not
    :param cx_smiles_fields: flags for cxsmiles creation
    :param normalize_dative_bonds: converts potential dative bonds drawn as single bonds to dative bonds
    :param include_query_features: if query features should be part of the hash or not
    :param seed: seed for xxhash. must be positive integer
    :return: a unique hash of the rdkit molecule
    """

    # Logic:
    # - Separate components/fragments
    # - Take cxsmiles of canonical tautomer without coordinates. These include enhanced stereo!
    # - Add query feature description to cxsmiles
    # - Repeat for each fragment
    # - concat all fragments and stereo code
    # - generate hash
    #
    # Part 1
    # GetMolFrags prior to rdkit 2023.03.1 loses enhanced stereo if there is more than one component in the molecule.
    components = separate_components(mol)

    # Part 2
    # Get canonical tautomer
    # Take cxsmiles of canonical tautomer and add query features to it
    # First we remove all conformers as we don't want coordinates in the cxsmiles
    # we want canonical smiles therefore the query features must be re-mapped to canonical smiles
    # atom indexes
    logger.debug("Generating hash for all components...")
    component_hashes = []
    for component in components:

        # iterating the python list is faster than lookup by idx and much faster than GetAtoms()
        # hence we iterate once to get it into a list and reuse that list for future iterations
        has_query_atom = False
        has_query_bond = False
        atoms = []
        query_atoms = []

        for i in range(0, component.GetNumAtoms()):
            atom = component.GetAtomWithIdx(i)
            atoms.append(atom)
            if atom.HasQuery():
                has_query_atom = True
                # Fix Query atoms:
                # A query list like [F,Cl,Br] leads to a SMILES containing "F", the first atom in the
                # list instead of a "*" when this same query is read from SMARTS.
                query_description = atom.DescribeQuery()
                if "AtomOr" in query_description: # likely too simplistic
                    logger.info(f"Found AtomOr query with atomic number {atom.GetAtomicNum()}. Setting it to 0.")
                    atom.SetAtomicNum(0)
                # Fix atomic number with atypical query order
                # A SMARTS input of [R1&C] will set the atom in SMILES as wildcard instead of a Carbon
                # because it is set to atomic number 0. Fix this and set to the proper AtomType
                if atom.GetAtomicNum() == 0 and query_description.count("AtomType") == 1:
                    atomic_num = int(p_atom_type.search(query_description).group(1))
                    atom.SetAtomicNum(atomic_num)
                query_atoms.append(atom)

        bonds = []
        query_bonds = []
        for i in range(0, component.GetNumBonds()):
            bond = component.GetBondWithIdx(i)
            bonds.append(bond)
            if bond.HasQuery():
                has_query_bond = True
                # Set bond to specific type
                # a SMARTS of [!@;:] will lead to an unspecified bond type while [:;!@] (different order) results in
                # aromatic bond type. If only 1 specific bond type in query, set it to that bond type
                query_description = bond.DescribeQuery()
                if query_description.count("BondOrder") == 1:
                    bond_type = int(p_bond_type.search(query_description).group(1))
                    bond.SetBondType(Chem.BondType.values[bond_type])
                # either or query means bond is unspecified
                elif query_description.count("BondOrder") > 1:
                    bond.SetBondType(Chem.BondType.UNSPECIFIED)
                query_bonds.append(bond)

        if normalize_dative_bonds:
            component = _normalize_dative_bonds(component, atoms, bonds)

        if not tautomer_sensitive:
            tauts = enumerator.Enumerate(component)
            if len(tauts) > 1:
                logger.debug("Found more than 1 tautomer. Using canonical tautomer.")
                canon_mol = enumerator.Canonicalize(component)
            else:
                canon_mol = component
        else:
            canon_mol = component

        write_params = Chem.SmilesWriteParams() # default write params
        component_hash = Chem.MolToCXSmiles(canon_mol, params=write_params, flags=cx_smiles_fields)

        # remove any text data stored in sgroups (SgD:Text) which has no chemical meaning
        # no option to omit this specifically, other sgroup data might be chemically relevant
        component_hash = re.sub(r"SgD:Text:.+::::", "", component_hash)
        # removing text might lead to empty additional data for cx smiles, remove that as well
        # Example: CC(=O)OC1CCCCC1 |,SgD:Text:XYZ::::,SgD:Text:ABC::::| => CC(=O)OC1CCCCC1 |,| => CC(=O)OC1CCCCC1
        component_hash = re.sub(r"\|[, ]?\|", "", component_hash).strip()

        if include_query_features and (has_query_atom or has_query_bond):
            # if molecule contains query atoms or bonds, special handling is needed to ensure constant hash
            # it requires that atoms need to be reordered in the order of the smiles output

            # _smilesAtomOutputOrder property is needed to map query bonds correctly
            # a copy of the molecule is created to match the atom indexes of the output order
            # for some reason GetPropsAsDict takes a significant amount of time (200µs) therefore
            # it's worth it to have this check. it also prevents iteration over bonds an atoms.
            # in total this if block reduces runtime by about 40% !!! (v2022.09.05)

            # GetPropsAsDict has unpredictable performance and can be very slow if the molecule contains many large
            # properties. For example loading from cxsmiles with coordinates increases the call by about 130µs.
            # For predictable performance we get the property directly as string and convert to list with
            # ast.literal_eval which is safe (read the docs)
            # order = canon_mol.GetPropsAsDict(True, True)["_smilesAtomOutputOrder"]
            o = canon_mol.GetProp("_smilesAtomOutputOrder")
            order = ast.literal_eval(o)
            logger.debug(order)
            # working on a reordered copy guarantees the same order for same molecules
            #m2 = Chem.RenumberAtoms(canon_mol, order)
            # newOrder is [3,2,0,1], then atom 3 in the original molecule will be atom 0 in the new one
            # re order the list so it works for sorting zipped lists
            sort_order = [None] * len(order)
            for idx, val in enumerate(order):
                sort_order[val] = idx

            # working on a reordered atom list guarantees the same order for same molecules
            # this is much faster than Chem.RenumberAtoms(canon_mol, order), 5 vs 40 µs
            # see:
            # https://stackoverflow.com/questions/6618515/sorting-list-according-to-corresponding-values-from-a-parallel-list
            # https://www.pythoncentral.io/how-to-sort-a-list-tuple-or-object-with-sorted-in-python/
            # create list of lists, sort list by first element, the new order, return the second element into new list
            atoms = [x for _, x in sorted(zip(sort_order, atoms), key=lambda item: item[0])]

            logger.debug("Hashing Query features...")
            bonds_hashed = []
            # for atom in m2.GetAtoms():
            # https://github.com/rdkit/rdkit/issues/6208 GetAtoms() is slow
            for atom in atoms:
                logger.debug(f"Atom Index: {atom.GetIdx()}")
                for bond in atom.GetBonds():
                    if bond.HasQuery() and bond.GetIdx() not in bonds_hashed:
                        q = _canonicalize_query(bond.GetSmarts())
                        b = bond.GetBeginAtomIdx()
                        e = bond.GetEndAtomIdx()
                        # map to canonical atom order
                        b = order.index(b)
                        e = order.index(e)
                        # always use the lowest atom index first
                        srt = sorted([b,e])
                        component_hash += ' |{}:{},{}|'.format(q, srt[0], srt[1])
                        bonds_hashed.append(bond.GetIdx())

                if atom.HasQuery():
                    qry = _canonicalize_atom_query(atom.GetSmarts())
                    atom_idx = order.index(atom.GetIdx())
                    component_hash += ' |{}:{}|'.format(atom_idx, qry)

        component_hashes.append(component_hash)

    # canonical component order
    component_hashes.sort()
    h = xxhash.xxh3_64(seed=seed)
    for ch in component_hashes:
        logger.debug(f"Raw hash for component: {ch}.")
        h.update(ch.encode('ASCII'))
    hex_hash = h.hexdigest()
    logger.debug(f"Generated hash {hex_hash} for molecule {Chem.MolToSmiles(mol)}.")
    return hex_hash


def _canonicalize_query(smarts: str):
    add_braces = False
    if smarts.startswith("["):
        smarts = smarts.strip("[]")
        add_braces = True
    can_smarts = ";".join(sorted(
        ",".join(sorted(
            "&".join(sorted(factor.split("&")))
            for factor in term.split(",")
        ))
        for term in smarts.split(";")
    ))
    if add_braces:
        can_smarts = "[" + can_smarts + "]"
    return can_smarts


def _canonicalize_atom_query(smarts: str):
    return re.sub(r'(\[[a-zA-Z0-9#,;:&]+\])', lambda m: _canonicalize_query(m.group()), smarts)


def _normalize_dative_bonds(mol: Chem.RWMol, atoms: list, bonds: list) -> Chem.RWMol:

    mol = _single_to_dative_bonds(mol, atoms)
    mol = _remove_dative_bonds(mol, bonds)
    return mol


def _is_transition_metal(atom: Chem.Atom) -> bool:
    n = atom.GetAtomicNum()
    return (22 <= n <= 29) or (40 <= n <= 47) or (72 <= n <= 79)


def _single_to_dative_bonds(mol: Chem.RWMol, atoms: list, from_atoms=(7, 8)) -> Chem.RWMol:
    """
    Replaces single bonds between metals and atoms with atomic numbers in fom_atoms
    with dative bonds. The replacement is only done if the atom has "too many" bonds.

    see https://www.rdkit.org/docs/Cookbook.html (source of this code)

    :param mol: molecule to replace single bonds with dative bonds
    :param from_atoms: source atomic numbers of the single bonds to consider replacement
    :return: the modified molecule with dative bonds

    """
    pt = Chem.GetPeriodicTable()
    mol.UpdatePropertyCache(strict=False)
    # https://github.com/rdkit/rdkit/issues/6208 GetAtoms() is slow hence lookup by index
    metals = [atom for atom in atoms if _is_transition_metal(atom)]
    for metal in metals:
        for nbr in metal.GetNeighbors():
            if nbr.GetAtomicNum() in from_atoms and \
                    nbr.GetExplicitValence() > pt.GetDefaultValence(nbr.GetAtomicNum()) and \
                    mol.GetBondBetweenAtoms(nbr.GetIdx(), metal.GetIdx()).GetBondType() == Chem.BondType.SINGLE:
                mol.RemoveBond(nbr.GetIdx(), metal.GetIdx())
                mol.AddBond(nbr.GetIdx(), metal.GetIdx(), Chem.BondType.DATIVE)
    return mol


def _remove_dative_bonds(mol: Chem.RWMol, bonds: list) -> Chem.Mol:

    to_remove = []
    # https://github.com/rdkit/rdkit/issues/6208 GetAtoms() is slow and same applies to GetBonds()
    for bond in bonds:
        if bond.GetBondType() == Chem.BondType.DATIVE:
            begin_atom_idx = bond.GetBeginAtomIdx()
            end_atom_idx = bond.GetEndAtomIdx()
            to_remove.append((begin_atom_idx, end_atom_idx, bond))
    for dative_bond in to_remove:
        mol.RemoveBond(dative_bond[0], dative_bond[1])
        bonds.remove(dative_bond[2])
    # reassign double bond stereochemistry from coordinates
    if len(mol.GetConformers()) > 0:
        Chem.SetDoubleBondNeighborDirections(mol, mol.GetConformer(0))
    return mol
