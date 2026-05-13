#!/usr/bin/env python3

from os.path import dirname, join as pjoin
from yaml import safe_load
import unittest
from mcnp_utilities.lib.materials import Isotope, Element, Material, mix_materials, zaid_to_ZA, get_isotope_mass


with open(pjoin(dirname(__file__), '..', 'data', 'iso_ref.yaml'), 'r') as f:
  iref = safe_load(f)


class TestIsotope(unittest.TestCase):
  def test_initialization(self):
    isotope = Isotope('8016')
    self.assertEqual(isotope.name, 'Oxygen-16')
    self.assertEqual(isotope.Z, 8)
    self.assertEqual(isotope.A, 16)

    isotope = Isotope('H-1')
    self.assertEqual(isotope.name, 'Hydrogen-1')
    self.assertEqual(isotope.Z, 1)
    self.assertEqual(isotope.A, 1)

    isotope = Isotope('C-12')
    self.assertEqual(isotope.name, 'Carbon-12')
    self.assertEqual(isotope.Z, 6)
    self.assertEqual(isotope.A, 12)

    with self.assertRaises(AssertionError):
      isotope = Isotope('8000')

    with self.assertRaises(AssertionError):
      isotope = Isotope('346000')

    with self.assertRaises(AssertionError):
      isotope = Isotope('Q-47')

    with self.assertRaises(AssertionError):
      isotope = Isotope('Quarkion-47')

  def test_mass(self):
    isotope = Isotope('8016')
    self.assertAlmostEqual(isotope.get_mass(), iref[8]['isotopes'][16]['atomic mass'], places=4)

  def test_name(self):
    isotope = Isotope('8016')
    self.assertEqual(str(isotope), 'Oxygen-16')


class TestElement(unittest.TestCase):
  def test_initialization(self):
    element = Element(8)
    self.assertEqual(element.Z, 8)
    for A in iref[8]['isotopes']:
      self.assertIn(A, element.isotope_list)
      self.assertAlmostEqual(element.atom_fractions[A], iref[8]['isotopes'][A]['atom fraction'])
      self.assertAlmostEqual(element.weight_fractions[A], iref[8]['isotopes'][A]['mass fraction'])
    self.assertEqual(element.name, iref[8]['name'])
    self.assertEqual(element.zaid, '8000')
    self.assertFalse(element.exclusion_list)

    example_fs = {16 : 0.5, 17 : 0.5}
    with self.assertRaises(AssertionError):
      element = Element(8, atom_fractions=example_fs, weight_fractions=example_fs)

    with self.assertRaises(AssertionError):
      element = Element(435)

    element = Element(8, atom_fractions=example_fs)
    self.assertAlmostEqual(element.get_mass(), 0.5 * iref[8]['isotopes'][16]['atomic mass'] + 0.5 * iref[8]['isotopes'][17]['atomic mass'])

    element = Element(8, weight_fractions=example_fs)
    self.assertAlmostEqual(element.get_mass(), 1.0 / sum([f / iref[8]['isotopes'][A]['atomic mass'] for A, f in example_fs.items()]))

  def test_mass(self):
    element = Element(8)
    self.assertAlmostEqual(element.get_mass(), 15.999, places=3)

  def test_name(self):
    element = Element(8)
    self.assertEqual(str(element), 'Oxygen')


class TestMaterial(unittest.TestCase):
  def test_initialization_list(self):
    nucs = {'6012' : 0.5, '8016' : 0.5}
    material = Material('atom', nuclides=nucs)
    self.assertEqual(material.number, 1)
    for k in nucs.keys():
      self.assertAlmostEqual(material.nuclide_atom_fractions[k], 0.5)
    for k, v in nucs.items():
      self.assertAlmostEqual(material.nuclide_weight_fractions[k], v * get_isotope_mass(k) / material.get_molar_mass())
    self.assertAlmostEqual(material.get_molar_mass(), 0.5 * iref[6]['isotopes'][12]['atomic mass'] + 0.5 * iref[8]['isotopes'][16]['atomic mass'])

  def test_initialization_add(self):
    nucs = {'6012' : 0.5, '8016' : 0.5}
    material = Material('atom')
    for k, v in nucs.items():
      material.add_nuclide(Isotope(k), v)
    material.normalize()
    self.assertEqual(material.number, 1)
    for k in nucs.keys():
      self.assertAlmostEqual(material.nuclide_atom_fractions[k], 0.5)
    for k, v in nucs.items():
      self.assertAlmostEqual(material.nuclide_weight_fractions[k], v * get_isotope_mass(k) / material.get_molar_mass())
    self.assertAlmostEqual(material.get_molar_mass(), 0.5 * iref[6]['isotopes'][12]['atomic mass'] + 0.5 * iref[8]['isotopes'][16]['atomic mass'])

  def test_split(self):
    material = Material('atom', nuclides={'8000' : 1.0})
    material.split_elements()
    self.assertEqual(len(material.nuclides), len(iref[8]['isotopes']))
    for k, v in material.nuclide_atom_fractions.items():
      Z, A = zaid_to_ZA(k)
      self.assertAlmostEqual(v, iref[Z]['isotopes'][A]['atom fraction'])
    for k, v in material.nuclide_weight_fractions.items():
      Z, A = zaid_to_ZA(k)
      self.assertAlmostEqual(v, iref[Z]['isotopes'][A]['mass fraction'])

  def test_combine(self):
    nucs = {'8016' : 0.5, '8017' : 0.5}
    material = Material('atom', nuclides=nucs)
    material.combine_isotopes()
    self.assertEqual(len(material.nuclides), 1)
    self.assertEqual(material.nuclides[0].zaid, '8000')
    self.assertAlmostEqual(material.get_molar_mass(), sum([v * get_isotope_mass(k) for k, v in nucs.items()]))

  def test_mix(self):
    m1 = Material('atom', nuclides={'6012' : 0.45, '6013' : 0.55})
    m2 = Material('atom', nuclides={'8016' : 0.15, '8017' : 0.85})
    new_mat = mix_materials([m1, m2], [0.1, 0.9], 'atom', num=2)
    for zaid, f in zip(['6012', '6013', '8016', '8017'], [0.45 * 0.1, 0.55 * 0.1, 0.15 * 0.9, 0.85 * 0.9]):
      self.assertAlmostEqual(new_mat.nuclide_atom_fractions[zaid], f)
    self.assertEqual(new_mat.number, 2)

if __name__ == '__main__':
  unittest.main()
