#!/usr/bin/env python3


import unittest
import numpy as np
from mcnp_utilities.utils.paramesh import cell


class TestParamesh(unittest.TestCase):

  test_points = np.array([
    [0, 0, 0],
    [1, 0, 0],
    [0, 1, 0],
    [1, 1, 0],
    [0, 0, 1],
    [1, 0, 1],
    [0, 1, 1],
    [1, 1, 1]
  ])
  test_cyl = np.array([
    [1, 0, 0],
    [2, 0, 0],
    [1, 0, 0.25],
    [2, 0, 0.25],
    [1, 1, 0],
    [2, 1, 0],
    [1, 1, 0.25],
    [2, 1, 0.25]
  ])
  test_val = 10.5

  def test_cell_init(self):
    test_cell = cell(self.test_points, self.test_val, self.test_val)
    self.assertIsNotNone(test_cell.points)
    self.assertIsNotNone(test_cell.values)

  def test_convert_to_Cyl(self):
    test_cell = cell(self.test_cyl, self.test_val, self.test_val)

    test_cell.convert_to_rec()

    answers = np.array([
      [1, 0, 0],
      [2, 0, 0],
      [0, 1, 0],
      [0, 2, 0],
      [1, 0, 1],
      [2, 0, 1],
      [0, 1, 1],
      [0, 2, 1]
    ])
    for i, pt in enumerate(answers):
      for j, val in enumerate(pt):
        self.assertAlmostEqual(val, answers[i, j])


if __name__ == '__main__':
  unittest.main()
