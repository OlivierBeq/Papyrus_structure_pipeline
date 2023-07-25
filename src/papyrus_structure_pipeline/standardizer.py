# -*- coding: utf-8 -*-

#
#  Copyright (c) 2023 Olivier J. M. Béquignon
#  All rights reserved.
#
#  This file is part of the Papyrus_StructurePipeline project.
#  The contents are covered by the terms of the MIT license
#  which is included in the file LICENSE, found at the root
#  of the source tree.

"""Standardizing methods."""

from enum import Enum, auto
from typing import Optional, Tuple, Union

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize.rdMolStandardize import TautomerEnumerator, Uncharger
from rdkit.rdBase import BlockLogs

# Remove warnings of importing the ChEMBL Structure Pipeline
with BlockLogs():
    from chembl_structure_pipeline import standardizer
    from chembl_structure_pipeline.exclude_flag import exclude_flag


class StandardizationResult(Enum):
    """Result of the standardization process."""
    CORRECT_MOLECULE = auto()
    NON_SMALL_MOLECULE = auto()
    INORGANIC_MOLECULE = auto()
    MIXTURE_MOLECULE = auto()
    STANDARDIZATION_ERROR = auto()


class InorganicSubtype(Enum):
    """Reason the compound is deemed inorganic."""
    IS_ORGANIC = auto()
    NO_CC_BOND = auto()
    NOT_CHONPSFClIBrB = auto()
    EXCLUSION_FLAG_SET = auto()  # molecule was flagged by the ChEMBL Structure Pipeline


class SaltStrippingResult(Enum):
    """Result of removing salts, metals and additional substructures."""
    STRIPPED_MOLECULE = auto()
    EMPTY_MOLECULE = auto()


# Define salts not flagged by the ChEMBL structure pipeline
SALTS = ['[Na+]', '[K+]', 'Cl', 'C1=NC=CC=C1', '[O-][Cl+3]([O-])([O-])[O-]', 'O=CN(C)C', '[No]',
         'C', 'C1COCCO1', 'CCCCCCN(C)C', 'CCCCCC', 'CN(C)C(=N)N', 'CN(C)C',
         'c1ccc([B-](c2ccccc2)(c2ccccc2)c2ccccc2)cc1']

# Define metals as all but H, B, C, N, O, P, F, S, Cl, Br, I
METALS = ["He", "Li", "Be", "Ne", "Na", "Mg", "Al", "Si", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co",
          "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
          "Ag", "Cd", "In", "Sn", "Sb", "Te", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
          "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi",
          "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md",
          "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"]
METALS = [f'[{metal}]' for metal in METALS]

# Define organic atoms
ORGANIC_ATOMS = ['C', 'H', 'O', 'N', 'P', 'S', 'F', 'Cl', 'I', 'Br']


def standardize(mol: Chem.Mol,
                remove_additional_salts: bool = True,
                remove_additional_metals: bool = True,
                filter_mixtures: bool = True,
                filter_inorganic: bool = True,
                uncharge: bool = True,
                filter_non_small_molecule: bool = True,
                canonicalize_tautomer: bool = True,
                small_molecule_min_mw: float = 200,
                small_molecule_max_mw: float = 800,
                tautomer_allow_stereo_removal: bool = True,
                tautomer_max_tautomers: int = 2 ** 32 - 1,
                return_type: bool = False,
                raise_error: bool = True
                ) -> Union[Chem.Mol, Tuple[Optional[Chem.Mol], StandardizationResult]]:
    """Standardize a molecule from either a RDKit molecule or SMILES.

    Steps:
    1) Use the ChEMBL structure pipeline (get parent molecule & standardization)
    2) remove supplementary salts & metals
    3) check if a mixture
    4) check if organic
    5) uncharge
    5) check if 200 <= MW <= 800
    6) canonicalize tautomer
    7) Use the ChEMBL structure pipeline (get parent molecule & standardization)

    :param mol: RDKit molecule to be standardized
    :param remove_additional_salts: should additional salts then those
    dealt with in the ChEMBL Structure Pipeline be removed
    :param remove_additional_metals: should salt removal also remove metal atoms
    :param filter_mixtures: should mixtures be filtered out
    :param filter_inorganic: should inorganic molecules be filtered out
    :param uncharge: should formal charges be removed (based on a list of substructures)
    :param filter_non_small_molecule: should non-small molecules be filtered out
    :param canonicalize_tautomer: should a canonical tautomer be determined
    :param small_molecule_min_mw: minimum molecular weight to be a valid small molecule
    :param small_molecule_max_mw: maximum molecular weight to be a valid small molecule
    :param tautomer_allow_stereo_removal: allow stereocenters to be removed by tautomerization
    :param tautomer_max_tautomers: maximum number of tautomers to enumerate (ignored if 0)
    :param return_type: If True, include the `StandardizationResult` in the return value
    :param raise_error: If True, standardization failures raise exceptions, if False, None is returned
    (with the standardization error if `return_type` is set).
    :return: a tuple of either the standardized molecule or None and a standardization_result flag
    """
    # Apply ChEMBL standardization
    try:
        std_mol = _apply_chembl_standardization(mol)
        # Remove additional salts
        if remove_additional_salts:
            std_mol = _remove_supplementary_salts(std_mol, include_metals=remove_additional_metals)
        # Ensure is not a mixture
        if filter_mixtures and is_mixture(std_mol):
            return (None, StandardizationResult.MIXTURE_MOLECULE) if return_type else None
        # Ensure organic
        if filter_inorganic and not is_organic(std_mol):
            return (None, StandardizationResult.INORGANIC_MOLECULE) if return_type else None
        # Remove charges
        if uncharge:
            std_mol = _uncharge(std_mol)
        # Ensure small molecule
        if filter_non_small_molecule and not is_small_molecule(std_mol,
                                                               min_molwt=small_molecule_min_mw,
                                                               max_molwt=small_molecule_max_mw):
            return (None, StandardizationResult.NON_SMALL_MOLECULE) if return_type else None
        # Obtain canonical tautomer
        if canonicalize_tautomer:
            std_mol = _canonicalize_tautomer(std_mol,
                                             allow_stereo_removal=tautomer_allow_stereo_removal,
                                             max_tautomers=tautomer_max_tautomers)
        # Apply standardization once again
        final_mol = _apply_chembl_standardization(std_mol)
        return (final_mol, StandardizationResult.CORRECT_MOLECULE) if return_type else final_mol
    except Exception as e:
        if raise_error:
            raise RuntimeError('Molecule could not be standardized') from e
        elif return_type:
            return (None, StandardizationResult.STANDARDIZATION_ERROR)
        else:
            return None


def _apply_chembl_standardization(mol: Chem.Mol) -> Chem.Mol:
    """Apply the ChEMBL structure standardization pipeline to a RDKit molecule.

    :param mol: RDKit molecule to be standardized
    """
    if mol is None:
        raise ValueError('Either RDKit molecule or SMILES must be specified')
    # Use the ChEMBL structure pipeline to standardize the molecule
    block = BlockLogs()  # Disable RDKit outputs
    standardized_mol, _ = standardizer.get_parent_mol(standardizer.standardize_mol(mol))
    standardized_smiles = Chem.MolToSmiles(standardized_mol)
    del block  # Re-enable them
    # Ensure standardized SMILES can be read back
    try:
        block = BlockLogs()  # Disable RDKit outputs
        standardized_mol = Chem.MolFromSmiles(standardized_smiles)
        del block  # Re-enable them
        if standardized_mol is None:
            raise ValueError(f'Could not parse standardized SMILES: {standardized_smiles}')
    except Exception as e:
        raise ValueError(f'Could not parse standardized SMILES: {standardized_smiles}') from e
    return standardized_mol


def _canonicalize_tautomer(mol: Chem.Mol, allow_stereo_removal: bool = True,
                           max_tautomers: int = 2 ** 32 - 1) -> Chem.Mol:
    """Obtain the RDKit canonical tautomer of the given RDKit molecule."""
    if mol is None:
        raise ValueError('A RDKit molecule must be specified')
    block = BlockLogs()  # Disable RDKit outputs
    # Parameter of tautomer enumeration
    enumerator = TautomerEnumerator()
    enumerator.SetMaxTautomers(max_tautomers)
    # Enumerate tautomers
    if allow_stereo_removal:
        enumerator.SetRemoveSp3Stereo(True)
        enumerator.SetRemoveBondStereo(True)
        tautos = enumerator.Enumerate(mol)
    else:
        enumerator.SetRemoveSp3Stereo(False)
        enumerator.SetRemoveBondStereo(False)
        tautos = enumerator.Enumerate(mol)
        # Keep tautomers with stereo unchanged if required
        def get_num_ciral_centers(mol_: Chem.Mol):
            """Utility function."""
            return len([atom.GetChiralTag()
                        for atom in mol_.GetAtoms()
                        if atom.GetChiralTag() != AllChem.ChiralType.CHI_UNSPECIFIED])

        orig_chiral_centers = get_num_ciral_centers(mol)
        if orig_chiral_centers > 0:
            tautos = [tauto for tauto in tautos
                      if get_num_ciral_centers(tauto) == orig_chiral_centers]
    # Determine canonical tautomer
    can_tauto = enumerator.PickCanonical(tautos)
    del block  # Re-enable them
    if can_tauto is not None:
        return can_tauto
    raise ValueError(f'Could not obtain canonical tautomer: {Chem.MolToSmiles(mol)}')


def is_organic(mol: Chem.Mol, return_type: bool = False) -> Union[bool, Tuple[bool, InorganicSubtype]]:
    """Return whether the RKDit molecule is organic or not.

    Makes use of the `ORGANIC_ATOMS` variable to identify inorganic atoms.

    :param mol: RDKit molecule
    :param return_type: If True, include the `InorganicSubtype` in the return value
    :return: True if return_type is True and the molecule has no exclusion flag
    set by the `chembl_structure_pipeline`, has at least one C-C bond and is only
    made of C, H, O, N, P, S and halogen atoms.
    """
    if mol is None:
        raise ValueError('A RDKit molecule must be specified')
    # Obtain ChEMBL structure exclusion flag
    if exclude_flag(Chem.MolToMolBlock(mol)):
        return False if not return_type else (False, InorganicSubtype.EXCLUSION_FLAG_SET)
    # Ensure at least 1 C-C bond
    ccbond = Chem.MolFromSmarts('[#6]~[#6]')
    if len(mol.GetSubstructMatches(ccbond)) == 0:
        return False if not return_type else (False, InorganicSubtype.NO_CC_BOND)
    # Ensure contains only C, H, O, N, P, S and halogens
    global ORGANIC_ATOMS
    if len(set(atom.GetSymbol() for atom in mol.GetAtoms()).difference(set(ORGANIC_ATOMS))):
        return False if not return_type else (False, InorganicSubtype.NOT_CHONPSFClIBrB)
    return True if not return_type else (True, InorganicSubtype.IS_ORGANIC)


def is_small_molecule(mol: Chem.Mol,
                      min_molwt: float = 200,
                      max_molwt: float = 800,
                      ) -> bool:
    """Return whether the RDKit molecule is a 'small molecule'.

    :param mol: molecule to examine
    :param min_molwt: molecular weight under which the molecule is not considered a small molecule
    :param max_molwt: molecular weight above which the molecule is not considered a small molecule
    """
    if mol is None:
        raise ValueError('A RDKit molecule must be specified')
    return min_molwt <= AllChem.CalcExactMolWt(mol) <= max_molwt


def is_mixture(mol: Chem.Mol) -> bool:
    """Return whether the molecule is composed of multiple fragments."""
    if mol is None:
        raise ValueError('A RDKit molecule must be specified')
    return len(Chem.GetMolFrags(mol)) > 1


def _remove_supplementary_salts(mol: Chem.Mol,
                                include_metals: bool = True,
                                return_type: bool = False
                                ) -> Union[Chem.Mol, Tuple[Chem.Mol, SaltStrippingResult]]:
    """Remove salts not dealt with by the ChEMBL pipeline structure.

    Adapted from the ChEMBL structure pipeline.
    https://github.com/chembl/ChEMBL_Structure_Pipeline.git

    Makes use of the `SALTS` and `METALS` variables to identify substructures to remove.

    :param include_metals: should the variable `METALS` be used toremove additional substructures.
    :param return_type: If True, include the `SaltStrippingResult` in the return value
    :return: the stripped molecule
    """
    if mol is None:
        raise ValueError('A RDKit molecule must be specified')
    # Parse salts
    salts = (SALTS + METALS) if include_metals else SALTS
    for i, _ in enumerate(salts):
        salts[i] = Chem.MolFromSmiles(salts[i])
    # Obtain fragments as molecules
    frags = list(Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False))
    # Sanitize fragments
    for i, _ in enumerate(frags):
        frag = Chem.RemoveHs(frags[i], sanitize=False)
        frag.UpdatePropertyCache(strict=False)
        Chem.SetAromaticity(frag)
        frags[i] = frag
    # Flag fragments to keep
    keep = [1] * len(frags)
    for i, frag in enumerate(frags):
        for salt in salts:
            if (keep[i] and frag.GetNumAtoms() == salt.GetNumAtoms()
                    and frag.GetNumBonds() == salt.GetNumBonds() and frag.HasSubstructMatch(salt)
                    and salt.HasSubstructMatch(frag)):
                keep[i] = 0
            if not max(keep):
                # All was removed, return initial molecule
                return Chem.Mol(mol) if not return_type else (Chem.Mol(mol), SaltStrippingResult.EMPTY_MOLECULE)
    # Collect fragments to keep
    frags = [frags[i] for i, x in enumerate(keep) if x == 1]
    # Only one fragment left
    if len(frags) == 1:
        return frags[0] if not return_type else (frags[0], SaltStrippingResult.STRIPPED_MOLECULE)
    # Otherwise combine remaining fragments
    else:
        res = frags[0]
        for frag in frags[1:]:
            res = Chem.CombineMols(res, frag)
        return res if not return_type else (res, SaltStrippingResult.STRIPPED_MOLECULE)


def _uncharge(mol: Chem.Mol) -> Chem.Mol:
    """Uncharge the given RDKit molecule.

    Adapted from Hans de Winter's contribution to B3DB
    (https://github.com/theochem/B3DB/blob/main/cleaning/02_clean_smiles_chembl_way_20210215.py)
    """
    if mol is None:
        raise ValueError('A RDKit molecule must be specified')
    uncharger = Uncharger()
    mol = uncharger.uncharge(mol)
    # Do a second pass on zwitterions
    romol = Chem.RWMol(mol)
    neutralizing_smiles = [  # Imidazoles
        ('[n+;H]', 'n'),
        # Amines
        ('[N+;!H0]', 'N'),
        # Carboxylic acids and alcohols
        ('[$([O-]);!$([O-][#7])]', 'O'),
        # Thiols
        ('[S-;X1]', 'S'),
        # Sulfonamides
        ('[$([N-;X2]S(=O)=O)]', 'N'),
        # Enamines
        ('[$([N-;X2][C,N]=C)]', 'N'),
        # Tetrazoles
        ('[n-]', '[nH]'),
        # Sulfoxides
        ('[$([S-]=O)]', 'S'),
        # Amides
        ('[$([N-]C=O)]', 'N')
    ]
    neutralizing_reactions = [(Chem.MolFromSmarts(x), Chem.MolFromSmiles(y, False)) for x, y in neutralizing_smiles]
    for i, (reactant, product) in enumerate(neutralizing_reactions):
        while romol.HasSubstructMatch(reactant):
            rms = AllChem.ReplaceSubstructs(romol, reactant, product)
            romol = rms[0]
    mol = Chem.Mol(romol)
    Chem.SanitizeMol(mol)
    return mol
