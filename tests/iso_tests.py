#!/usr/bin/env python3

from argparse import ArgumentTypeError
from mcnp_utilities.lib.materials import Material
from mcnp_utilities.utils.iso import determine_components, process_input, nonnegative_int
from unittest import TestCase, main

class TestIso(TestCase):
  def test_determine_components(self):
    """
    Docstring for test_determine_components

    Tests functionality of determine_components function.
    """
    comps = determine_components('H')
    self.assertEqual(len(comps), 1)
    self.assertEqual(comps[0], 'H')

    comps = determine_components('H,O')
    self.assertEqual(len(comps), 2)
    self.assertEqual(comps[0], 'H')
    self.assertEqual(comps[1], 'O')

    comps = determine_components('(H),O')
    self.assertEqual(len(comps), 2)
    self.assertEqual(comps[0], '(H)')
    self.assertEqual(comps[1], 'O')

    comps = determine_components('(H),(O)')
    self.assertEqual(len(comps), 2)
    self.assertEqual(comps[0], '(H)')
    self.assertEqual(comps[1], '(O)')

    comps = determine_components('H:1,O:2')
    self.assertEqual(len(comps), 2)
    self.assertEqual(comps[0], 'H:1')
    self.assertEqual(comps[1], 'O:2')

  def test_process_input(self):
    """
    Docstring for test_process_input

    Tests functionality of process_input function.
    """
    mats = process_input('H', 0, quiet=True)
    self.assertTrue(isinstance(mats, Material))
    self.assertEqual(mats.nuclides[0].zaid, '1000')
    self.assertAlmostEqual(mats.nuclide_atom_fractions['1000'], 1)

    mats = process_input('H:2,O:1', 0, quiet=True)
    self.assertTrue(isinstance(mats, Material))
    self.assertEqual(len(mats.nuclides), 2)
    self.assertEqual(mats.nuclides[0].zaid, '1000')
    self.assertEqual(mats.nuclides[1].zaid, '8000')
    self.assertAlmostEqual(mats.nuclide_atom_fractions['1000'], 2/3)
    self.assertAlmostEqual(mats.nuclide_atom_fractions['8000'], 1/3)

    mats = process_input('(H:1,O:1):0.5,(B:1,C:1):0.5', 0, quiet=True)
    self.assertTrue(isinstance(mats, Material))
    self.assertEqual(len(mats.nuclides), 4)
    self.assertAlmostEqual(mats.nuclide_atom_fractions['1000'], 0.25)
    self.assertAlmostEqual(mats.nuclide_atom_fractions['8000'], 0.25)
    self.assertAlmostEqual(mats.nuclide_atom_fractions['5000'], 0.25)
    self.assertAlmostEqual(mats.nuclide_atom_fractions['6000'], 0.25)

  def test_nonnegative_int(self):
    """
    Docstring for test_nonnegative_int

    Tests functionality of nonnegative_int function.
    """
    with self.assertRaises(ArgumentTypeError):
      nonnegative_int(-1)
    with self.assertRaises(ArgumentTypeError):
      nonnegative_int('test')

if __name__ == '__main__':
  main()
