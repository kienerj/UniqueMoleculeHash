from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
import copy
import xxhash
import rdkit


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
    components = Chem.GetMolFrags(mol, asMols=True, frags=frags, fragsMolAtomMapping=mapping)

    # GetMolFrags prior to 2023.03.1 loses enhanced stereo if there is more than one component (structure) in the
    # molecule. Therefore the enhanced stereo needs to be recreated and reassigned.

    if rdkit.__version__ >= "2023.03.01":
        return components

    # tuple to list for later assignments
    components = list(components)

    if len(components) > 1 and len(mol.GetStereoGroups()) > 0:

        stereo_groups = [[] for i in range(len(components))]

        for st_grp in mol.GetStereoGroups():
            atoms = st_grp.GetAtoms()
            # group is limited to 1 component
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

    return components


def get_unique_hash(mol: Chem.Mol) -> str:
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
    enumerator = rdMolStandardize.TautomerEnumerator()
    component_hashes = []
    for component in components:
        tauts = enumerator.Enumerate(component)
        if len(tauts) > 1:
            canon_mol = enumerator.Canonicalize(component)
        else:
            canon_mol = component

        # Sets _smilesAtomOutputOrder property which we need to map query bonds correctly
        Chem.MolToSmiles(canon_mol)
        # create copy and remove conformers so that cxsmiles doesn't contain coordinates!
        mol_copy = copy.deepcopy(canon_mol)
        mol_copy.RemoveAllConformers()
        component_hash = Chem.MolToCXSmiles(mol_copy)
        order = canon_mol.GetPropsAsDict(True, True)["_smilesAtomOutputOrder"]
        m2 = Chem.RenumberAtoms(mol_copy, order)
        bonds = []
        for atom in m2.GetAtoms():
            for bond in atom.GetBonds():
                if bond.HasQuery() and bond.GetIdx() not in bonds:
                    bonds.append(bond.GetIdx())
                    q = bond.GetSmarts()
                    b = bond.GetBeginAtomIdx()
                    e = bond.GetEndAtomIdx()
                    component_hash += ' |{}:{},{}|'.format(q, b, e)
        component_hash += ';'
        component_hashes.append(component_hash)

    # canonical component order
    component_hashes.sort()
    h = xxhash.xxh64()
    for ch in component_hashes:
        h.update(ch.encode('ASCII'))
    return h.hexdigest()