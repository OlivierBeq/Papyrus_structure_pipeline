
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Papyrus Structure Pipeline

Papyrus protocols used to standardize molecules. First used in Papyrus 05.6
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7377161.svg)](https://doi.org/10.5281/zenodo.7377161).

## Installation

From source:

    git clone https://github.com/OlivierBeq/Papyrus_Structure_Pipeline.git
    pip install ./Papyrus_Structure_Pipeline

with pip:

```bash
pip install papyrus_structure_pipeline
```

## Usage

### Standardize a compound

Comparison to the [ChEMBL Structure Pipeline](https://github.com/chembl/ChEMBL_Structure_Pipeline):

```python
from rdkit import Chem
from chembl_structure_pipeline import standardizer as ChEMBL_standardizer
from papyrus_structure_pipeline import standardizer as Papyrus_standardizer

# CHEMBL1560279
smiles = "CCN(CC)C(=O)[n+]1ccc(OC)cc1.c1ccc([B-](c2ccccc2)(c2ccccc2)c2ccccc2)cc1"

mol = Chem.MolFromSmiles(smiles)
out1 = ChEMBL_standardizer.standardize_mol(mol)
out2 = Papyrus_standardizer.standardize(mol)


print(Chem.MolToSmiles(out1))
# CCN(CC)C(=O)[n+]1ccc(OC)cc1.c1ccc([B-](c2ccccc2)(c2ccccc2)c2ccccc2)cc1

print(Chem.MolToSmiles(out2))
# CCN(CC)C(=O)[n+]1ccc(OC)cc1
```

Get details on the standardization to identify why it fails for some molecules:

```python
smiles_list = [
    # erlotinib
    "n1cnc(c2cc(c(cc12)OCCOC)OCCOC)Nc1cc(ccc1)C#C",
    # midecamycin
    "CCC(=O)O[C@@H]1CC(=O)O[C@@H](C/C=C/C=C/[C@@H]([C@@H](C[C@@H]([C@@H]([C@H]1OC)O[C@H]2[C@@H]([C@H]([C@@H]([C@H](O2)C)O[C@H]3C[C@@]([C@H]([C@@H](O3)C)OC(=O)CC)(C)O)N(C)C)O)CC=O)C)O)C",
    # selenofolate
    "C1=CC(=CC=C1C(=O)NC(CCC(=O)OCC[Se]C#N)C(=O)O)NCC2=CN=C3C(=N2)C(=O)NC(=N3)N",
    # cisplatin
    "N.N.Cl[Pt]Cl"
]

for smiles in smiles_list:
    mol = Chem.MolFromSmiles(smiles)
    print(Papyrus_standardizer.standardize(mol, return_type=True))

    
# (<rdkit.Chem.rdchem.Mol object at 0x000000946F99B580>, <StandardizationResult.CORRECT_MOLECULE: 1>)
# (None, <StandardizationResult.NON_SMALL_MOLECULE: 2>)
# (None, <StandardizationResult.INORGANIC_MOLECULE: 3>)
# (None, <StandardizationResult.MIXTURE_MOLECULE: 4>)
```

Allow other atoms to be considered organic:

```python
smiles = "CCN(CC)C(=O)C1=CC=C(S1)C2=C3C=CC(=[N+](C)C)C=C3[Se]C4=C2C=CC(=C4)N(C)C.F[P-](F)(F)(F)(F)F"
mol = Chem.MolFromSmiles(smiles)

print(Papyrus_standardizer.standardize(mol, return_type=True))
# (None, <StandardizationResult.INORGANIC_MOLECULE: 3>)

Papyrus_standardizer.ORGANIC_ATOMS.append('Se')

print(Papyrus_standardizer.standardize(mol, return_type=True))
# (<rdkit.Chem.rdchem.Mol object at 0x0000009F24D15F90>, <StandardizationResult.CORRECT_MOLECULE: 1>)

Papyrus_standardizer.ORGANIC_ATOMS = Papyrus_standardizer.ORGANIC_ATOMS[:-1]

print(Papyrus_standardizer.standardize(mol, return_type=True))
# (None, <StandardizationResult.INORGANIC_MOLECULE: 3>)
```

Add custom substructures to be removed as salts:

```python
# lomitapide
smiles = "C1CN(CCC1NC(=O)C2=CC=CC=C2C3=CC=C(C=C3)C(F)(F)F)CCCCC4(C5=CC=CC=C5C6=CC=CC=C64)C(=O)NCC(F)(F)F.c1ccccc1"
mol = Chem.MolFromSmiles(smiles)

print(Papyrus_standardizer.standardize(mol, return_type=True))
# (None, <StandardizationResult.MIXTURE_MOLECULE: 4>)

Papyrus_standardizer.SALTS.append('c1ccccc1')

print(Papyrus_standardizer.standardize(mol, return_type=True))
# (<rdkit.Chem.rdchem.Mol object at 0x0000009F24D15F90>, <StandardizationResult.CORRECT_MOLECULE: 1>)

Papyrus_standardizer.SALTS = Papyrus_standardizer.SALTS[:-1]

print(Papyrus_standardizer.standardize(mol, return_type=True))
# (None, <StandardizationResult.MIXTURE_MOLECULE: 4>)
```


## Documentation


```python
def standardize(mol,
                remove_additional_salts=True, remove_additional_metals=True,
                filter_mixtures=True, filter_inorganic=True, filter_non_small_molecule=True,
                canonicalize_tautomer=True, small_molecule_min_mw=200, small_molecule_max_mw=800,
                tautomer_allow_stereo_removal=True, tautomer_max_tautomers=0, return_type=False
                ) -> Chem.Mol:
```

#### Parameters

- ***mol  : Chem.Mol***  
    RDKit molecule object to standardize.
- ***remove_additional_salts  : bool***  
    Removes a custom set of fragments if present in the molecule object.
- ***remove_additional_metals  : bool***  
    Removes metal fragments if present in the molecule object.<br/>Ignored if `remove_additional_salts` is set to `False`.
- ***filter_mixtures  : bool***  
    Return `None` if the molecule is a mixture.  
- ***filter_inorganic  : bool***  
    Return `None` if the molecule is a inorganic.
- ***filter_non_small_molecule  : bool***  
    Return `None` if the molecule is not a small molecule.
- ***canonicalize_tautomer  : bool***  
    Canonicalize the tautomeric state of the molecule.
- ***small_molecule_min_mw  : float***  
    Molecular weight under which a molecule is considered too small. 
- ***small_molecule_max_mw  : float***  
    Molecular weight above which a molecule is considered too big.
- ***tautomer_allow_stereo_removal  : bool***  
    Allow the tautomer search algorithm to remove stereocenters. 
- ***tautomer_max_tautomers  : int (< 2 <sup>32</sup>)***  
    Maximum number of tautomers to consider by the tautomer search algorithm.
- ***return_type  : bool***  
    Add a StandardizationResult to the return value. 

________________

```python
def is_organic(mol, return_type=False) -> bool:
```

Return whether the RKDit molecule is organic or not.
 - Makes internal use of the variable `ORGANIC_ATOMS` 

#### Parameters

- ***mol  : Chem.Mol***  
    RDKit molecule object to check the organic nature of.
- ***return_type  : bool***  
    Add a InorganicSubtype to the return value. 

________________

```python
def is_small_molecule(mol,
                      min_molwt=200,
                      max_molwt=800
                      ) -> bool:
```

Return whether the RKDit molecule has a molecular weight within `min_molwt` and `max_molwt`.

#### Parameters

- ***mol  : Chem.Mol***  
    RDKit molecule object to check the organic nature of.
- ***min_molwt  : float***  
    Molecular weight under which a molecule is considered too small. 
- ***max_molwt  : float***  
    Molecular weight above which a molecule is considered too big.

________________


```python
def is_mixture(mol) -> bool:
```

Return whether the RKDit molecule is composed of multiple fragments.

#### Parameters

- ***mol  : Chem.Mol***  
    RDKit molecule object to check the organic nature of.

________________

```python
def _apply_chembl_standardization(mol) -> Chem.Mol:
``` 

Apply the ChEMBL structure standardization pipeline on a RDKit molecule.
- Makes use of both `ChEMBL_standardizer.get_parent_mol` and `ChEMBL_standardizer.standardize_mol`.
- Ensures the obtained SMILES can be parsed by the RDKit.

#### Parameters

- ***mol  : Chem.Mol***  
    RDKit molecule object to apply the ChEMBL Structure Pipeline to.

________________

```python
def _canonicalize_tautomer(mol,
                           allow_stereo_removal=True, max_tautomers=2 ** 32 - 1
                           ) -> Chem.Mol:
``` 

Obtain the RDKit canonical tautomer of the given RDKit molecule.
- Makes use of both `ChEMBL_standardizer.get_parent_mol` and `ChEMBL_standardizer.standardize_mol`.
- Ensures the obtained SMILES can be parsed by the RDKit.

#### Parameters

- ***mol  : Chem.Mol***  
    RDKit molecule object to RDKit canonical tautomer of.
- ***allow_stereo_removal  : bool***  
    Allow the tautomer search algorithm to remove stereocenters.
- ***max_tautomers  : int (<2 <sup>32</sup>)***  
    Maximum number of tautomers to consider by the tautomer search algorithm.

________________

```python
def _remove_supplementary_salts(mol,
                                include_metals=True, return_type=False
                                ) -> Chem.Mol:
``` 

Remove substructures (e.g. salts) not dealt with by the ChEMBL pipeline structure.

The additional substructures are defined by the `SALTS` variable.

#### Parameters

- ***mol  : Chem.Mol***  
    RDKit molecule object from which to remove additional substructures.
- ***include_metals  : bool***  
    Removes metal fragments if present in the molecule object (defined by the `METALS` variable).
- ***return_type  : bool***  
    Add a SaltStrippingResult to the return value.

