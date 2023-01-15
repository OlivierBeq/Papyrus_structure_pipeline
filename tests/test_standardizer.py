# -*- coding: utf-8 -*-

"""Tests for Descriptor."""


import sys
import unittest

from papyrus_structure_pipeline import standardizer
from tests.constants import *
from rdkit.Chem.AllChem import ChiralType


class TestStandardizer(unittest.TestCase):
    """Tests for standardizer."""

    def test_standardizer_correct(self):
        self.assertIsNotNone(standardizer.standardize(MOLECULES['CHEMBL1560279']))
        self.assertIsNotNone(standardizer.standardize(MOLECULES['erlotinib']))
        self.assertIsNotNone(standardizer.standardize(MOLECULES['lomitapide']))

    def test_standardizer_non_small(self):
        self.assertEqual(standardizer.standardize(MOLECULES['midecamycin'], return_type=True)[1],
                         standardizer.StandardizationResult.NON_SMALL_MOLECULE)

    def test_standardizer_inorganic(self):
        self.assertEqual(standardizer.standardize(MOLECULES['selenofolate'], return_type=True)[1],
                         standardizer.StandardizationResult.INORGANIC_MOLECULE)
        self.assertEqual(standardizer.standardize(MOLECULES['CHEMBL457061'], return_type=True)[1],
                         standardizer.StandardizationResult.INORGANIC_MOLECULE)

    def test_standardizer_mixture(self):
        # standardizer.SALTS.extend(['c1ccccc1', ])
        self.assertEqual(standardizer.standardize(MOLECULES['cisplatin'], return_type=True)[1],
                         standardizer.StandardizationResult.MIXTURE_MOLECULE)

    def test_standardizer_non_small_pass(self):
        self.assertIsNotNone(standardizer.standardize(MOLECULES['midecamycin'], filter_non_small_molecule=False))
        self.assertIsNotNone(standardizer.standardize(MOLECULES['midecamycin'], small_molecule_max_mw=820))

    def test_standardizer_inorganic_pass(self):
        self.assertIsNotNone(standardizer.standardize(MOLECULES['selenofolate'], filter_inorganic=False))
        self.assertIsNotNone(standardizer.standardize(MOLECULES['CHEMBL457061'], filter_inorganic=False))
        standardizer.ORGANIC_ATOMS.append('Se')
        self.assertIsNotNone(standardizer.standardize(MOLECULES['selenofolate']))
        self.assertIsNotNone(standardizer.standardize(MOLECULES['CHEMBL457061']))
        standardizer.ORGANIC_ATOMS = standardizer.ORGANIC_ATOMS[:-1]
        self.assertEqual(standardizer.standardize(MOLECULES['selenofolate'], return_type=True)[1],
                         standardizer.StandardizationResult.INORGANIC_MOLECULE)
        self.assertEqual(standardizer.standardize(MOLECULES['CHEMBL457061'], return_type=True)[1],
                         standardizer.StandardizationResult.INORGANIC_MOLECULE)

    def test_standardizer_mixture_pass(self):
        self.assertEqual(standardizer.standardize(MOLECULES['cisplatin'], filter_mixtures=False, return_type=True)[1],
                         standardizer.StandardizationResult.INORGANIC_MOLECULE)
        self.assertIsNotNone(standardizer.standardize(MOLECULES['cisplatin'], filter_mixtures=False, filter_inorganic=False))

    def test_standardizer_tautomer(self):
        mol = Chem.MolFromSmiles('OC1=NC(=O)[C@H](C=C1)N1C(=O)C2=CC=CCC2=C1O')
        std_mol = standardizer.standardize(mol, canonicalize_tautomer=False)
        self.assertEqual(Chem.MolToSmiles(mol), Chem.MolToSmiles(std_mol))

        std_mol = standardizer.standardize(mol, canonicalize_tautomer=True)
        self.assertNotEqual(Chem.MolToSmiles(mol), Chem.MolToSmiles(std_mol))

    def test_standardizer_tautomer_stereocenter(self):
        mol = Chem.MolFromSmiles('OC1=NC(=O)[C@H](C=C1)N1C(=O)C2=CC=CCC2=C1O')
        std_mol = standardizer.standardize(mol, filter_non_small_molecule=False, tautomer_allow_stereo_removal=True)
        self.assertFalse('@' in Chem.MolToSmiles(std_mol))

        std_mol = standardizer.standardize(mol, filter_non_small_molecule=False, tautomer_allow_stereo_removal=False)
        self.assertTrue('@' in Chem.MolToSmiles(std_mol))
