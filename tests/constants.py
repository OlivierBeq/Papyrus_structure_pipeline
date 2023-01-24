# -*- coding: utf-8 -*-

"""Constants for unit tests."""

from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover


MOLECULES = {
    'CHEMBL1560279': Chem.MolFromSmiles("CCN(CC)C(=O)[n+]1ccc(OC)cc1"),
    'erlotinib': Chem.MolFromSmiles("n1cnc(c2cc(c(cc12)OCCOC)OCCOC)Nc1cc(ccc1)C#C"),
    'midecamycin': Chem.MolFromSmiles("CCC(=O)O[C@@H]1CC(=O)O[C@@H](C/C=C/C=C/[C@@H]([C@@H](C[C@@H]([C@@H]([C@H]1OC)O"
                                      "[C@H]2[C@@H]([C@H]([C@@H]([C@H](O2)C)O[C@H]3C[C@@]([C@H]([C@@H](O3)C)OC(=O)CC)"
                                      "(C)O)N(C)C)O)CC=O)C)O)C"),
    'selenofolate': Chem.MolFromSmiles("C1=CC(=CC=C1C(=O)NC(CCC(=O)OCC[Se]C#N)C(=O)O)NCC2=CN=C3C(=N2)C(=O)NC(=N3)N"),
    'cisplatin': Chem.MolFromSmiles("N.N.Cl[Pt]Cl"),
    'CHEMBL457061': Chem.MolFromSmiles("CCN(CC)C(=O)C1=CC=C(S1)C2=C3C=CC(=[N+](C)C)C=C3[Se]C4=C2C=CC(=C4)N(C)C"),
    'lomitapide': Chem.MolFromSmiles("C1CN(CCC1NC(=O)C2=CC=CC=C2C3=CC=C(C=C3)C(F)(F)F)CCCCC4(C5=CC=CC=C5C6=CC=CC=C64)"
                                     "C(=O)NCC(F)(F)F")
}

ADDITIONAL_SALTS = [Chem.MolFromSmiles("c1ccc([B-](c2ccccc2)(c2ccccc2)c2ccccc2)cc1"),
                    Chem.MolFromSmiles("F[P-](F)(F)(F)(F)F"),
                    Chem.MolFromSmiles("c1ccccc1")] + SaltRemover().salts
