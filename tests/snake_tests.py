#!/usr/bin/env python3

from os import chdir, getcwd, mkdir
from os.path import dirname, exists, join as pjoin
from shutil import rmtree
from mcnp_utilities.utils.snake import snake, SNAKE_DEFAULTS
from unittest import TestCase, main


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

class TestSnake(TestCase):
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
    for fyle, expected_content in test_answers.items():
      self.assertTrue(exists(pjoin(getcwd(), fyle)))
      with open(fyle, 'r') as f:
        lines = f.readlines()
        if not isinstance(expected_content, list):
          raise TypeError("Expected content must be a list")
        self.assertEqual(len(expected_content), len(lines))
        for i, expected_line in enumerate(expected_content):
          self.assertEqual(expected_line.strip(), lines[i].strip())

  def test_external_usage(self):
    """
    Docstring for test_external_usage

    Tests that an external library can be used.
    """
    args = InputArgs(
      pjoin(input_dir, 'test-external-usage.snake'),
      lib=pjoin(input_dir, 'test-external-usage-lib.py')
    )
    snake(args, quiet=True)

    self.checkAnswers({
      'test-external-usage_0_0' : ['1 a const tsnoc'],
      'test-external-usage_0_1' : ['1 b const tsnoc'],
      'test-external-usage_1_0' : ['2 a const tsnoc'],
      'test-external-usage_1_1' : ['2 b const tsnoc']
    })

  def test_resolve_error(self):
    """
    Docstring for test_resolve_error

    Ensures an error is thrown when a variable
    cannot be resolved.
    """
    args = InputArgs(pjoin(input_dir, 'test-resolve-error.snake'))
    with self.assertRaises(RecursionError):
      snake(args, quiet=True)

  def test_external_library(self):
    """
    Docstring for test_external_library

    Tests that external library call works and that
    variables can be assigned to variables that are
    returned from an external function call.
    """
    args = InputArgs(
      pjoin(input_dir, 'test-external-library.snake'),
      lib=pjoin(input_dir, 'test-external-library-lib.py')
    )
    snake(args, quiet=True)

    self.checkAnswers({
      'test-external-library_0' : ['1 -1 -1'],
      'test-external-library_1' : ['2 -2 -2']
    })

  def test_comment_character(self):
    """
    Docstring for test_comment_character

    Tests different comment character from default
    and f-string printing.
    """
    args = InputArgs(
      pjoin(input_dir, 'test-comment-character.snake'),
      cchar='!'
    )
    snake(args, quiet=True)

    self.checkAnswers({
      'test-comment-character_0' : ['x = 0.250'],
      'test-comment-character_1' : ['x = 0.500'],
      'test-comment-character_2' : ['x = 0.750'],
      'test-comment-character_3' : ['x = 1.000']
    })

  def test_multi_assignment(self):
    """
    Docstring for test_multi_assignment

    Tests multi-variable assignment.
    """
    args = InputArgs(pjoin(input_dir, 'test-multi-assignment.snake'))
    snake(args, quiet=True)

    self.checkAnswers({
      'test-multi-assignment_0' : ['10 20 30']
    })

  def test_external_library_keywords(self):
    """
    Docstring for test_external_library_keywords

    Tests that keyword arguments in external library functions are handled correctly.
    """
    args = InputArgs(
      pjoin(input_dir, 'test-external-library-keywords.snake'),
      lib=pjoin(input_dir, 'test-external-library-keywords-lib.py')
    )
    snake(args, quiet=True)

    self.checkAnswers({
      'test-external-library-keywords_0' : ['(1, 14, 6)'],
      'test-external-library-keywords_1' : ['(2, 14, 6)'],
      'test-external-library-keywords_2' : ['(3, 14, 6)'],
    })

  def test_f_strings(self):
    """
    Docstring for test_f_strings

    Tests that f-strings with nested braces are handled correctly.
    """
    args = InputArgs(pjoin(input_dir, 'test-f-strings.snake'))
    snake(args, quiet=True)

    self.checkAnswers({
      'test-f-strings_0' : ['the value of x is: 1, this is a string containing x: 1'],
      'test-f-strings_1' : ['the value of x is: 2, this is a string containing x: 2'],
      'test-f-strings_2' : ['the value of x is: 3, this is a string containing x: 3'],
    })

  def test_multi_pass(self):
    """
    Docstring for test_multi_pass

    Tests constants that are correctly evaluated on a second pass.
    """
    args = InputArgs(pjoin(input_dir, 'test-multi-pass.snake'))
    snake(args, quiet=True)

    self.checkAnswers({
      'test-multi-pass_0' : ['value = 3']
    })

  def test_noniterable_key(self):
    """
    Docstring for test_noniterable_key

    Tests that non-iterable key values raise a TypeError.
    """
    args = InputArgs(pjoin(input_dir, 'test-noniterable-key.snake'))

    with self.assertRaises(TypeError):
      snake(args, quiet=True)

  def test_material_mixing(self):
    """
    Docstring for test_material_mixing

    Tests manual material specification and mixing.
    """
    args = InputArgs(pjoin(input_dir, 'test-material-mixing.snake'))
    snake(args, quiet=True)

    self.checkAnswers({
      'test-material-mixing_0' : [
        'C Material 5',
        'M5    8000  1.000000E+00 $ Oxygen'
      ],
      'test-material-mixing_1' : [
        'C Material 5',
        'M5    1000  5.000000E-01 $ Hydrogen',
        '      8000  5.000000E-01 $ Oxygen'
      ],
      'test-material-mixing_2' : [
        'C Material 5',
        'M5    1000  1.000000E+00 $ Hydrogen'
      ]
    })

  def test_compendium_mixing(self):
    """
    Docstring for test_compendium_mixing

    Tests material generation from PNNL compendium and mixing.
    """
    args = InputArgs(pjoin(input_dir, 'test-compendium-mixing.snake'))
    snake(args, quiet=True)

    self.checkAnswers({
      'test-compendium-mixing_0' : [
        'C Material 1',
        'M1    1001  -5.593641E-02 $ Hydrogen-1',
        '      1002  -1.285697E-05 $ Hydrogen-2',
        '      8016  -5.019535E-01 $ Oxygen-16',
        '      8017  -2.032116E-04 $ Oxygen-17',
        '      8018  -1.160765E-03 $ Oxygen-18',
        '     92234  -1.170396E-04 $ Uranium-234',
        '     92235  -1.322207E-02 $ Uranium-235',
        '     92236  -6.033259E-05 $ Uranium-236',
        '     92238  -4.273338E-01 $ Uranium-238'
      ],
    })

  def test_keyword_variable(self):
    """
    Docstring for test_keyword_variable

    Tests that a keyword variable name raises a ValueError.
    """
    args = InputArgs(pjoin(input_dir, 'test-keyword-variable.snake'))

    with self.assertRaises(ValueError):
      snake(args, quiet=True)

if __name__ == '__main__':
  main()
