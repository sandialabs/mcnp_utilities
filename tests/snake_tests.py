#!/usr/bin/env python3

from os import chdir, getcwd, mkdir
from os.path import dirname, exists, join as pjoin
from shutil import rmtree
from mcnp_utilities.utils.snake import snake, SNAKE_DEFAULTS
import unittest


input_dir = pjoin(dirname(__file__), 'snake_inputs')
snake_test_dir = pjoin(dirname(__file__), 'snake_testing')

class InputArgs:
  def __init__(self, fname, lib=None, cchar=SNAKE_DEFAULTS['COMMENT_CHARACTER']):
    self.input = fname
    self.library = lib
    self.keyfile = None
    self.print_only = SNAKE_DEFAULTS['PRINT_ONLY']
    self.override = True
    self.naming = SNAKE_DEFAULTS['NAMING_CONVENTION']
    self.delimiter = SNAKE_DEFAULTS['DELIMITER']
    self.comment_char = cchar
    self.extension = SNAKE_DEFAULTS['FILE_EXTENSION']
    self.organize = SNAKE_DEFAULTS['ORGANIZE']

class TestSnake(unittest.TestCase):
  def setUp(self) -> None:
    """
    Docstring for setUp

    Creates the snake working test directory
    and changes the working directory to it.
    """
    if not exists(snake_test_dir):
      mkdir(snake_test_dir)
    chdir(snake_test_dir)

  def tearDown(self) -> None:
    """
    Docstring for tearDown

    Changes working directory back and deletes
    snake working test directory.
    """
    chdir(dirname(__file__))
    rmtree(snake_test_dir)

  def checkAnswers(self, test_answers):
    """
    Docstring for checkAnswers

    Checks that correct files exist and contain the correct content.

    :param test_answers: dict containing files as keys and contents as answers
    """
    for fyle in test_answers.keys():
      self.assertTrue(exists(pjoin(getcwd(), fyle)))
      with open(fyle, 'r') as f:
        lines = f.readlines()
        self.assertEqual(len(lines), 1)
        self.assertEqual(test_answers[fyle], lines[0].strip())

  def test_basic_1(self):
    args = InputArgs(
      pjoin(input_dir, 'test-1.snake'),
      lib=pjoin(input_dir, 'test-1-lib.py')
    )
    snake(args, quiet=True)

    self.checkAnswers({
      'test-1_0_0' : '1 a const tsnoc',
      'test-1_0_1' : '1 b const tsnoc',
      'test-1_1_0' : '2 a const tsnoc',
      'test-1_1_1' : '2 b const tsnoc'
    })

  def test_basic_2(self):
    """
    Docstring for test_basic_2

    Ensures an error is thrown when a variable
    cannot be resolved.
    """
    args = InputArgs(pjoin(input_dir, 'test-2.snake'))
    with self.assertRaises(RecursionError):
      snake(args, quiet=True)

  def test_basic_3(self):
    """
    Docstring for test_basic_3

    Tests that external library call works and that
    variables can be assigned to variables that are
    returned from an external function call.
    """
    args = InputArgs(
      pjoin(input_dir, 'test-3.snake'),
      lib=pjoin(input_dir, 'test-3-lib.py')
    )
    snake(args, quiet=True)

    self.checkAnswers({
      'test-3_0' : '1 -1 -1',
      'test-3_1' : '2 -2 -2'
    })

  def test_basic_4(self):
    """
    Docstring for test_basic_4

    Tests different comment character from default
    and f-string printing.
    """
    args = InputArgs(
      pjoin(input_dir, 'test-4.snake'),
      cchar='!'
    )
    snake(args, quiet=True)

    self.checkAnswers({
      'test-4_0' : 'x = 0.250',
      'test-4_1' : 'x = 0.500',
      'test-4_2' : 'x = 0.750',
      'test-4_3' : 'x = 1.000'
    })

  def test_basic_5(self):
    """
    Docstring for test_basic_5

    Tests multi-variable assignment.
    """
    args = InputArgs(pjoin(input_dir, 'test-5.snake'))
    snake(args, quiet=True)

    self.checkAnswers({'test-5_0' : '10 20 30'})

if __name__ == '__main__':
  unittest.main()
