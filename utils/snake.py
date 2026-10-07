#!/usr/bin/env python3

# Snake (Python implementation of Worm)

# Inspired by the Worm utility developed by
# Tom Jones (LA-CC-99-69). The goal was to
# create a modern tool that can perform the
# same functions as Worm.

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from copy import copy
from numpy import where, linspace, logspace, arange, asarray
from math import pi, e, log, log10, sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, prod
from itertools import product
from importlib import import_module
from sys import path
from os import linesep
from os.path import realpath, dirname, split as psplit, join as pjoin, splitext
from datetime import datetime
from random import random, randint
from re import compile
from collections.abc import Iterable
# Local modules
from mcnp_utilities.lib.materials import get_compendium_material, get_compendium_material_card, mix_materials, Material
from mcnp_utilities.lib.basic_tools import create_nested_path


SNAKE_DEFAULTS = {
  'DELIMITER'         : '_',
  'COMMENT_CHARACTER' : '#',
  'FILE_EXTENSION'    : '',
  'NAMING_CONVENTION' : 'i',
  'OVERRIDE'          : False,
  'ORGANIZE'          : False,
  'PRINT_ONLY'        : False
}

constants = {
  'cm'     : 1,         'g'      : 1,             'mm'     : 0.1,           'kg'     : 1000,
  'm'      : 100,       'lb'     : 453.59237,     'inch'   : 2.54,          'oz'     : 28.349523125,
  'ft'     : 30.48,     'yd'     : 91.44,         'rad'    : 1,             'mil'    : 0.00254,
  'deg'    : pi / 180,  'cc'     : 1,             's'      : 1,             'l'      : 1000,
  'min'    : 60,        'ml'     : 1,             'hr'     : 3600,          'gal'    : 3785.411784,
  'day'    : 86400,     'ozfl'   : 29.5735295625, 'yr'     : 31556925.9747, 'bit'    : 0.0001,
  'e'      : e,         'an'     : 0.60221367,    'aN'     : 6.0221367E+23, 'pi'     : pi,
  'NOW'    : datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
  'DATE'   : datetime.now().strftime('%Y-%m-%d'),
  'TIME'   : datetime.now().strftime('%H:%M:%S')
}

allowable_fns = {
  'arange'   : arange, 'linspace' : linspace, 'logspace' : logspace,
  'abs'      : abs,    'log'      : log,      'log10'    : log10,
  'sin'      : sin,    'cos'      : cos,      'tan'      : tan,
  'asin'     : asin,   'acos'     : acos,     'atan'     : atan,
  'sinh'     : sinh,   'cosh'     : cosh,     'tanh'     : tanh,
  'random'   : random, 'randint'  : randint
}
material_fns = {
  'get_compendium_material'      : get_compendium_material,
  'get_compendium_material_card' : get_compendium_material_card,
  'mix_materials'                : mix_materials,
  'Material'                     : Material
}
allowable_fns.update(material_fns)

_ASSIGNMENT_LHS_RE = compile(
  r'^@?[A-Za-z_]\w*(?:\s*,\s*[A-Za-z_]\w*)*$'
)

def get_arguments():
  parser = ArgumentParser(description=r"""
  ______                    __
 /      \                  |  \
|  ▓▓▓▓▓▓\_______   ______ | ▓▓   __  ______
| ▓▓___\▓▓       \ |      \| ▓▓  /  \/      \
 \▓▓    \| ▓▓▓▓▓▓▓\ \▓▓▓▓▓▓\ ▓▓_/  ▓▓  ▓▓▓▓▓▓\
 _\▓▓▓▓▓▓\ ▓▓  | ▓▓/      ▓▓ ▓▓   ▓▓| ▓▓    ▓▓
|  \__| ▓▓ ▓▓  | ▓▓  ▓▓▓▓▓▓▓ ▓▓▓▓▓▓\| ▓▓▓▓▓▓▓▓
 \▓▓    ▓▓ ▓▓  | ▓▓\▓▓    ▓▓ ▓▓  \▓▓\\▓▓     \
  \▓▓▓▓▓▓ \▓▓   \▓▓ \▓▓▓▓▓▓▓\▓▓   \▓▓ \▓▓▓▓▓▓▓

--------------------------------------------
Scriptable Nesting and Keying Engine (SNAKE)
     (A Python implementation of WORM)
--------------------------------------------

This tool creates permutations of a file, with the
permutations determined by the cartesian product of key
variables. The permutations of the established key values
are used to create separate files in the cwd
(or specified directory structure).

Lines starting with # (by default) are not propagated to
generated files. These lines must be used to set the value
of constants and variables. Other lines are copied into the
generated files.

Constants and keys (variables) are set within curly braces.
Constants and keys can have any allowable Python variable
name (but keys must start with an @ symbol), and are set
using an equals sign. Constants may be scalars or lists
(created using Python list syntax). Keys must be lists.

Constants can be set using Python's multi-variable assignment syntax.

Examples:

{a, b = 1, 2}
{a, b = (1, 2)}
{a, b = [1, 2]}
{a, b = get_pair()}

The right-hand side must evaluate to an iterable with the same number
of values as the number of target variables.

This will not work for keys.

The following Python/NumPy functions can be used:
""" + f"{linesep}".join([f'- {k}' for k in allowable_fns if k not in material_fns]) + """

User-defined functions can also be imported.

The following material functions are available:
""" + f"{linesep}".join([f'- {k}(#)' for k in material_fns]) + """

The following built-in constants are available:
""" + f"{linesep}".join([f'- {k:4s} : {v}' for k, v in constants.items()]) + """

Values can be printed using Python f-string format codes.
A value is printed if there is no equals sign in the
curly braces and the line does not start with the comment
indicator character.

Files are named according to the index of the keys
(unless otherwise specified), in the order that the
keys are identified.

Example snake input file that will create 120 files:

# {@x = arange(1, 5, 1)}
# {@y = linspace(1, 20, 10)}
# {@z = [1, log(5/2), 10*2]}
# {c = [1, b, y]} {b=5}
c x={x} y={y:.3f} z={z:<.5f} c[1]={c[1]:^10g} b={b:>5n}""", formatter_class=RawDescriptionHelpFormatter)
  parser.add_argument(
    'input',
    type=str,
    help='Input file',
    metavar='<input file>'
  )
  parser.add_argument(
    '-c', dest='comment_char',
    default=SNAKE_DEFAULTS['COMMENT_CHARACTER'],
    help='Sets the comment indicator character (default: %(default)s)',
    type=str,
    metavar='<comment char.>'
  )
  parser.add_argument(
    '-l', dest='library',
    type=str,
    help='path to external python library',
    metavar='<path>'
  )
  parser.add_argument(
    '-p', dest='print_only',
    action='store_true',
    default=SNAKE_DEFAULTS['PRINT_ONLY'],
    help='Only print key values (without creating all files)'
  )
  parser.add_argument(
    '-y', dest='override',
    action='store_true',
    default=SNAKE_DEFAULTS['OVERRIDE'],
    help='Override confirmation prompt for creating files'
  )
  parser.add_argument(
    '-o', dest='organize',
    action='store_true',
    default=SNAKE_DEFAULTS['ORGANIZE'],
    help='Organize files into directories based on key index or name (based on -n option)'
  )
  parser.add_argument(
    '-k', dest='keyfile',
    type=str,
    help='Name of key file',
    metavar='<key file>'
  )
  parser.add_argument(
    '-d', dest='delimiter',
    type=str,
    default=SNAKE_DEFAULTS['DELIMITER'],
    help='Delimiter for file names (default: %(default)s)',
    metavar='<delimiter>'
  )
  parser.add_argument(
    '-e', dest='extension',
    type=str,
    help='Extension for produced files (default: %(default)s)',
    default=SNAKE_DEFAULTS['FILE_EXTENSION'],
    metavar='<extension>'
  )
  parser.add_argument(
    '-n', dest='naming',
    type=str,
    choices=['i', 'v'],
    default=SNAKE_DEFAULTS['NAMING_CONVENTION'],
    help='Naming convention for files. i=name by key index. v=name by key value. (default: %(default)s)',
    metavar='<i|v>'
  )
  return parser.parse_args()

def validate_permutations(permutations):
  for key, values in permutations.items():
    if isinstance(values, str) or not isinstance(values, Iterable):
      raise TypeError(
        f'Permutation key "{key}" must evaluate to a non-string iterable; got {values!r}'
      )

def evaluate(vname, string, vars, not_evaluated, new_allowable_fns):
  '''
  Evaluate an expression based on defined constants.
  If a variable is not defined, save it for processing on a second pass.
  '''
  try:
    result = eval(string, new_allowable_fns, vars)
  except NameError:
    not_evaluated[vname] = string
    return None

  if result is not None:
    return result

  not_evaluated[vname] = string
  return None

def find_top_level_assignment(expression):
  '''
  Return the index of a top-level assignment operator, or -1 if the
  expression is not an assignment.

  This ignores '=' inside function calls, lists, tuples, dicts, strings,
  and comparison operators like ==, <=, >=, !=.
  '''
  depth = 0
  quote = None
  escaped = False

  for i, char in enumerate(expression):
    if quote:
      if escaped:
        escaped = False
      elif char == '\\':
        escaped = True
      elif char == quote:
        quote = None
      continue

    if char in ('"', "'"):
      quote = char
      continue

    if char in '([{':
      depth += 1
      continue

    if char in ')]}':
      depth -= 1
      continue

    if char == '=' and depth == 0:
      prev_char = expression[i - 1] if i > 0 else ''
      next_char = expression[i + 1] if i + 1 < len(expression) else ''

      # Skip comparison-style operators.
      if prev_char in ('=', '!', '<', '>') or next_char == '=':
        continue

      return i

  return -1

def find_matching_brace(text, open_idx):
  '''
  Find the closing brace that matches text[open_idx].

  Braces inside quoted strings are ignored, so expressions like

    {mat_name = f'{x} wt. % water'}

  are parsed as one complete Snake expression instead of stopping at
  the inner f-string brace.
  '''
  quote = None
  escaped = False
  brace_depth = 0

  for i in range(open_idx + 1, len(text)):
    char = text[i]

    if quote:
      if escaped:
        escaped = False
      elif char == '\\':
        escaped = True
      elif char == quote:
        quote = None
      continue

    if char in ('"', "'"):
      quote = char
      continue

    if char == '{':
      brace_depth += 1
      continue

    if char == '}':
      if brace_depth == 0:
        return i
      brace_depth -= 1
      continue

  raise ValueError(f'Unmatched opening brace in line: {text.rstrip()}')

def process_line(
    text,
    permutations,
    deferred_constants,
    deferred_permutations,
    new_constants,
    new_allowable_fns
  ):
  '''
  Process a line if it contains expression(s).
  '''
  i = 0

  while i < len(text):
    if text[i] != '{':
      i += 1
      continue

    close_idx = find_matching_brace(text, i)
    expression = text[i + 1:close_idx]

    assignment_idx = find_top_level_assignment(expression)

    # Not an assignment; leave it alone so replace_line() can evaluate it later.
    if assignment_idx == -1:
      i = close_idx + 1
      continue

    name = expression[:assignment_idx].strip()
    expr = expression[assignment_idx + 1:].strip()

    # If the left side is not a valid Snake assignment target, skip it.
    if not _ASSIGNMENT_LHS_RE.match(name):
      i = close_idx + 1
      continue

    if name.startswith('@'):
      key_name = name[1:]

      permutations[key_name] = evaluate(
        key_name,
        expr,
        new_constants,
        deferred_permutations,
        new_allowable_fns
      )
    else:
      if ',' in name:
        names = [n.strip() for n in name.split(',')]
        values = evaluate(
          ', '.join(names),
          expr,
          new_constants,
          deferred_constants,
          new_allowable_fns
        )

        if values is not None:
          if len(values) != len(names):
            raise ValueError(
              f'Cannot unpack {len(values)} values into {len(names)} variables: {name} = {expr}'
            )

          for n, value in zip(names, values):
            new_constants[n] = value
        else:
          for j, n in enumerate(names):
            new_constants[n] = None
            deferred_constants[n] = f'({expr})[{j}]'
      else:
        new_constants[name] = evaluate(
          name,
          expr,
          new_constants,
          deferred_constants,
          new_allowable_fns
        )

    # Skip over the entire expression so we do not parse nested f-string
    # fields like {x} as separate Snake expressions.
    i = close_idx + 1

def replace_line(text, vars, new_allowable_fns):
  '''
  Replace instances of formatted variables in a line.
  '''
  while '{' in text:
    start_idx = text.index('{')
    close_idx = find_matching_brace(text, start_idx)
    v = text[start_idx:close_idx + 1]

    try:
      evaluated_variable = eval(f'f{v!r}', new_allowable_fns, vars)
    except NameError as err:
      var_name = v.strip('{}').split(':')[0]
      raise TypeError(f'Variable "{var_name}" is undefined!') from err
    except Exception as err:
      raise RuntimeError(f'Failed to evaluate output expression {v!r}') from err

    text = text[:start_idx] + str(evaluated_variable) + text[close_idx + 1:]

  return text

def write_keys(args, fobj, ps, n):
  '''
  Write the keys to the keyfile.
  '''
  fobj.write(f'{args.input}{2*linesep}')
  for k, v in ps.items():
    fobj.write(f'{k}  ({len(v)})  {",".join(str(val) for val in v)}{linesep}')
  fobj.write(f'{linesep}A total of {n} files.{2*linesep}')

def get_file_name(argo, indicies, kvs):
  '''
  Determine file name based on input arguments.
  '''
  ext = argo.extension
  if ext and not ext.startswith('.'):
    ext = '.' + ext

  base = splitext(psplit(argo.input)[-1])[0]

  if argo.naming == 'i':
    fname = f'{base}{argo.delimiter}{(argo.delimiter).join(indicies)}{ext if ext else ""}'
  elif argo.naming == 'v':
    fname = f'{base}{argo.delimiter}{(argo.delimiter).join((str(kv) for kv in kvs))}{ext if ext else ""}'

  if argo.organize:
    if argo.naming == 'i':
      key_path = pjoin(*indicies)
    elif argo.naming == 'v':
      key_path = pjoin(*(str(kv) for kv in kvs))
    fname = pjoin(key_path, fname)

  return fname

def write_file(args, fname, lines, file_number, total_files, vars, new_allowable_fns, quiet):
  '''
  Write a file with a specific set of key values.
  '''
  folder = dirname(fname)
  if folder:
    create_nested_path(folder)
  with open(fname, 'w') as f:
    for line in lines:
      if not line.startswith(args.comment_char):
        if '{' in line and '}' in line:
          f.write(replace_line(line, vars, new_allowable_fns))
        else:
          f.write(line)
    if not quiet:
      print(f'|>{(round(36*((file_number+1)/total_files))*"≈" + ">"):<37s}|', end='\r')

def snake(args, quiet=False):
  deferred_constants = {}
  deferred_permutations = {}
  permutations = {}
  key_file = None
  new_constants = copy(constants)
  new_allowable_fns = copy(allowable_fns)
  # Read external functions from specified file (if callable)
  if args.library:
    path.append(dirname(realpath(args.library)))
    new_allowable_fns.update({k : f for k, f in vars(import_module(psplit(args.library)[-1].rsplit('.')[0])).items() if callable(f) and not f.__name__.startswith("__")})
  # Read the lines in the file
  with open(args.input, 'r') as f:
    input_lines = f.readlines()
  # Read in all constants and lists of perumtation values
  [
    process_line(
      line,
      permutations,
      deferred_constants,
      deferred_permutations,
      new_constants,
      new_allowable_fns
    )
    for line in input_lines
    if '{' in line and '}' in line
  ]
  # Attempt to evaluate any expressions in keys that could not be evaluated on first pass (save expressions that are still unresolved, i.e., those depending on key values)
  resolved_any = True

  while resolved_any:
    resolved_any = False

    # First resolve deferred constants that only depend on other constants.
    for k, expr in list(deferred_constants.items()):
      try:
        result = eval(expr, new_allowable_fns, new_constants)
      except NameError:
        continue

      if result is not None:
        new_constants[k] = result
        del deferred_constants[k]
        resolved_any = True

    # Then resolve deferred permutation keys that depend on constants.
    for k, expr in list(deferred_permutations.items()):
      try:
        result = eval(expr, new_allowable_fns, new_constants)
      except NameError:
        continue

      if result is not None:
        permutations[k] = result
        del deferred_permutations[k]
        resolved_any = True

  if deferred_permutations:
    unresolved = ', '.join(
      f'{k} = {expr}' for k, expr in deferred_permutations.items()
    )
    raise NameError(
      f'The following key expressions could not be resolved before file generation: {unresolved}'
    )
  validate_permutations(permutations)

  # Print info and ask for confirmation (if not overridden)
  if not quiet:
    print(f'The following keys are identified:{linesep}')
    for k in permutations:
      print(f'- {k} ({len(permutations[k])} values)')
  nfiles = prod([len(v) for v in permutations.values()])
  if args.keyfile:
    key_file = open(args.keyfile, 'w')
    write_keys(args, key_file, permutations, nfiles)
  if args.print_only:
    if not quiet:
      print(f'{linesep}The key values are as follows:')
      for k in permutations:
        print(f'  {k}:')
        for v in permutations[k]:
          print(f'    - {v}')
    if args.keyfile:
      key_file.close()
    exit()
  if not quiet:
    print(f'{linesep}A total of {nfiles} files will be created.{linesep}')
  if not args.override:
    resp = input(f'Proceed [(y)/n]? ').lower()
    if resp and resp != 'y':
      print(f'{linesep}Stopping.')
      exit()
    else:
      print()
  indexed_values = [
    list(enumerate(values))
    for values in permutations.values()
  ]
  for fnum, combo in enumerate(product(*indexed_values)):
    idxs = [str(idx) for idx, value in combo]
    t = [value for idx, value in combo]
    # Copy existing constants
    iter_constants = copy(new_constants)
    # Add current permutation to variable dict
    for k, v in zip(permutations, t):
      iter_constants[k] = v
    # (Re)evaluate all expressions that depend on key values and update evaluator with new variables
    while any(v is None for v in iter_constants.values()):
      num_none = sum(v is None for v in iter_constants.values())
      for k, expr in deferred_constants.items():
        scratch_deferred = {}
        iter_constants[k] = evaluate(
          k,
          expr,
          iter_constants,
          scratch_deferred,
          new_allowable_fns
        )
      if sum(v is None for v in iter_constants.values()) == num_none:
        if not quiet:
          print('The following variables and expressions evaluate to None:')
          for name, none_expr in [
            (k, v)
            for k, v in deferred_constants.items()
            if evaluate(k, v, iter_constants, deferred_constants, new_allowable_fns) is None
          ]:
            print(f' - {name} : {none_expr}')
        raise RecursionError('Expressions are not being resolved! Ensure that user-defined functions cannot return None.')
    # Construct file name based on input arguments
    file_name = get_file_name(args, idxs, t)
    # Write current permutation to file
    write_file(args, file_name, input_lines, fnum, nfiles, iter_constants, new_allowable_fns, quiet)
    if args.keyfile:
      key_file.write(f'{file_name} {",".join(str(tv) for tv in t)}{linesep}')
  if args.keyfile:
    key_file.close()
  if not quiet:
    print(f'{linesep*2}All {nfiles} files successfully created.')

if __name__ == '__main__':
  snake(get_arguments())
