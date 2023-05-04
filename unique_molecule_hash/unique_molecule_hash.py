from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
import copy
import xxhash
import rdkit
import re
import logging


logger = logging.getLogger('unique_molecule_hash')


# Creating this instance is very costly so doing it only once instead of per function call reduces runtime by 5x!
tautomer_enumerator = rdMolStandardize.TautomerEnumerator()


def separate_components(mol: Chem.Mol) -> list:
    """
    Separates the rdkit mol into its individual components and returns them as a list of rdkit molecules.

    This function is required to circumvent a bug in rdkit before 2023.03.1. Chem.GetMolFrags removes enhancedstereo
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
    # molecule. Therefore the enhanced stereo needs to be recreated and reassigned.
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
                components[component_idx] = component.GetMol()
                logger.debug(f"Fixed Stereo for component {Chem.MolToSmiles(component)}")

    return components


# idea: make standard hash method and a configurable one to choose tautomer layer (inchi or enumerator), custom layer
# with or without enhanced stereo etc.
def get_unique_hash(mol: Chem.Mol, enumerator=tautomer_enumerator) -> str:
    """
    Creates a hash to compare RDKit molecules for uniqueness. it takes into account enhancedstereo, tautomerism and
    query features.

    Logic:
    - Separate components/fragments
    - Take cxsmiles of canonical tautomer without coordinates. These include enhanced stereo!
    - Add query feature description to cxsmiles
    - Repeat for each fragment
    - concat all fragments and stereo code
    - generate hash

    By default an internal rdMolStandardize.TautomerEnumerator() instance is used. You can pass your own as long
    as it has a "Canonicalize(mol)" method that returns a canonical tautomer. This of course impact the generated hash
    and also performance.

    :param mol: a valid rdkit molecule
    :param enumerator: tautomer enumerator to use
    :return: a unique hash of the rdkit molecule
    """

    # Part 1
    # GetMolFrags prior to rdkit 2023.03.1 loses enhanced stereo if there is more than one component in the molecule.
    components = separate_components(mol)

    # Part 2
    # Get canonical tautomer
    # Take cxmsiles of canonical tautomer and add query features to it
    # First we remove all conformers as we don't want coordinates in the cxsmiles
    # we want canonical smiles therefore the query features must be re-mapped to canonical smiles
    # atom indexes
    logger.debug("Generating hash for all components...")
    component_hashes = []
    for component in components:
        component = handle_dative_bonds(component)
        tauts = enumerator.Enumerate(component)
        if len(tauts) > 1:
            logger.debug("Found more than 1 tautomer. Using canonical tautomer.")
            canon_mol = enumerator.Canonicalize(component)
        else:
            canon_mol = component

        write_params = Chem.SmilesWriteParams() # default write params
        component_hash = Chem.MolToCXSmiles(canon_mol, params=write_params, flags=489)

        # remove any text data stored in sgroups (SgD:Text) which has no chemical meaning
        # no option to omit this specifically, other sgroup data might be chemically relevant
        component_hash = re.sub(r"SgD:Text:.+::::", "", component_hash)
        # removing text might lead to empty additional data for cx smiles, remove that as well
        # Example: CC(=O)OC1CCCCC1 |,SgD:Text:XYZ::::,SgD:Text:ABC::::| => CC(=O)OC1CCCCC1 |,| => CC(=O)OC1CCCCC1
        component_hash = re.sub(r"\|[, ]?\|", "", component_hash).strip()

        # if molecule contains query bonds, special handling is needed to ensure constant hash
        # it requires that atoms need to be reordered in the order of the smiles output
        if "~" in component_hash:
            # _smilesAtomOutputOrder property is needed to map query bonds correctly
            # a copy of the molecule is created to match the atom indexes of the output order
            # for some reason GetPropsAsDict takes a significant amount of time (200µs) therefore
            # it's worth it to have this check. it also prevents iteration over bonds an atoms.
            # in total this if block reduces runtime by about 40% !!! (v2022.09.05)
            order = canon_mol.GetPropsAsDict(True, True)["_smilesAtomOutputOrder"]
            m2 = Chem.RenumberAtoms(canon_mol, order)
            logger.debug("Determining Bond query features for hashing.")
            bonds = []
            for atom in m2.GetAtoms():
                for bond in atom.GetBonds():
                    if bond.HasQuery() and bond.GetIdx() not in bonds:
                        bonds.append(bond.GetIdx())
                        q = bond.GetSmarts()
                        b = bond.GetBeginAtomIdx()
                        e = bond.GetEndAtomIdx()
                        component_hash += ' |{}:{},{}|'.format(q, b, e)
        component_hashes.append(component_hash)

    # canonical component order
    component_hashes.sort()
    h = xxhash.xxh3_64()
    for ch in component_hashes:
        logger.debug(f"Raw hash for component: {ch}.")
        h.update(ch.encode('ASCII'))
    hex_hash = h.hexdigest()
    logger.debug(f"Generated hash {hex_hash} for molecule {Chem.MolToSmiles(mol)}.")
    return hex_hash


def handle_dative_bonds(mol: Chem.Mol) -> Chem.Mol:

    mol = single_to_dative_bonds(mol)
    mol = remove_dative_bonds(mol)
    return mol


def is_transition_metal(atom: Chem.Atom) -> bool:
    n = atom.GetAtomicNum()
    return (22 <= n <= 29) or (40 <= n <= 47) or (72 <= n <= 79)


def single_to_dative_bonds(mol: Chem.Mol, from_atoms=(7, 8)) -> Chem.Mol:
    """
    Replaces single bonds between metals and atoms with atomic numbers in fom_atoms
    with dative bonds. The replacement is only done if the atom has "too many" bonds.

    :param mol: molecule to replace single bonds with dative bonds
    :param from_atoms: source atomic numbers of the single bonds to consider replacement
    :return: the modified molecule with dative bonds

    """
    pt = Chem.GetPeriodicTable()
    rwmol = Chem.RWMol(mol)
    rwmol.UpdatePropertyCache(strict=False)
    metals = [at for at in rwmol.GetAtoms() if is_transition_metal(at)]
    for metal in metals:
        for nbr in metal.GetNeighbors():
            if nbr.GetAtomicNum() in from_atoms and \
                    nbr.GetExplicitValence() > pt.GetDefaultValence(nbr.GetAtomicNum()) and \
                    rwmol.GetBondBetweenAtoms(nbr.GetIdx(), metal.GetIdx()).GetBondType() == Chem.BondType.SINGLE:
                rwmol.RemoveBond(nbr.GetIdx(), metal.GetIdx())
                rwmol.AddBond(nbr.GetIdx(), metal.GetIdx(), Chem.BondType.DATIVE)
    return rwmol.GetMol()


def remove_dative_bonds(mol: Chem.Mol) -> Chem.Mol:

    rwmol = Chem.RWMol(mol)
    to_remove = []
    for bond in rwmol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DATIVE:
            begin_atom_idx = bond.GetBeginAtomIdx()
            end_atom_idx = bond.GetEndAtomIdx()
            to_remove.append((begin_atom_idx, end_atom_idx))
    for dative_bond in to_remove:
        rwmol.RemoveBond(dative_bond[0], dative_bond[1])

    # reassign double bond stereochemistry from coordiantes
    result =  rwmol.GetMol()
    if len(result.GetConformers()) > 0:
        Chem.SetDoubleBondNeighborDirections(result, result.GetConformer(0))
    return result
