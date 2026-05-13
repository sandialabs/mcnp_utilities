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
from os.path import realpath, dirname, split as psplit, join as pjoin
from datetime import datetime
from random import random, randint
from keyword import kwlist
# Local modules
from mcnp_utilities.lib.materials import material_fns
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
allowable_fns.update(**material_fns)

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

Constants can be set using Python's multi-variable
assignment syntax, but the right hand side MUST be
indexable. For example, the following assignment is valid:

{a, b = (1, 2)}

but the following is not:

{a, b = 1, 2}

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

def evaluate(vname, string, vars, not_evaluated, new_allowable_fns):
  '''
  Evaluate an expression based on defined constants.
  If a variable is not defined, save it for processing on a second pass.
  '''
  try:
    result = eval(string, new_allowable_fns, vars)
    if result is not None:
      return result
    else:
      not_evaluated[vname] = string
  except:
    not_evaluated[vname] = string

def process_line(text, permutations, not_evaluated, new_constants, new_allowable_fns):
  '''
  Process a line if it contains (an) expression(s).
  '''
  for i, char in enumerate(text):
    if char == '{':
      start_idx = i + 1
      end_idx = text[start_idx:].index('}')
      expression = text[start_idx:start_idx+end_idx]
      if '=' in expression:
        name, expr = (s.strip() for s in expression.split('='))
        if name.startswith('@'):
          permutations[name[1:]] = evaluate(name[1:], expr, new_constants, not_evaluated, new_allowable_fns)
        else:
          if ',' in name:
            for i, n in enumerate(name.split(',')):
              new_constants[n.strip()] = evaluate(n.strip(), f'{expr}[{i}]', new_constants, not_evaluated, new_allowable_fns)
          else:
            new_constants[name] = evaluate(name, expr, new_constants, not_evaluated, new_allowable_fns)

def replace_line(text, vars, not_evaluated, new_allowable_fns):
  '''
  Replace instances of formatted variables in a line.
  '''
  while '{' in text:
    start_idx = text.index('{')
    end_idx = text[start_idx:].index('}')+1
    v = text[start_idx:(start_idx+end_idx)]
    evaluated_variable = evaluate(v, f"f'{v}'", vars, not_evaluated, new_allowable_fns)
    if evaluated_variable is not None:
      text = text.replace(v, evaluated_variable)
    else:
      var_name = v.strip('{}').split(':')[0]
      if var_name in kwlist:
        raise ValueError(f'Variable "{var_name}" is a Python keyword and cannot be used as a variable name!')
      else:
        raise TypeError(f'Variable "{var_name}" is undefined!')
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
  if argo.naming == 'i':
    fname = f'{psplit(argo.input)[-1].rsplit(".")[0]}{argo.delimiter}{(argo.delimiter).join(indicies)}{"." + argo.extension if argo.extension else ""}'
  elif argo.naming == 'v':
    fname = f'{psplit(argo.input)[-1].rsplit(".")[0]}{argo.delimiter}{(argo.delimiter).join((str(kv) for kv in kvs))}{"." + argo.extension if argo.extension else ""}'
  if argo.organize:
    if argo.naming == 'i':
      key_path = pjoin(*indicies)
    elif argo.naming == 'v':
      key_path = pjoin(*(str(kv) for kv in kvs))
    fname = pjoin(key_path, fname)
  return fname

def write_file(args, fname, lines, file_number, total_files, vars, not_evaluated, new_allowable_fns, quiet):
  '''
  Write a file with a specific set of key values.
  '''
  create_nested_path(dirname(fname))
  with open(fname, 'w') as f:
    for line in lines:
      if not line.startswith(args.comment_char):
        f.write(replace_line(line, vars, not_evaluated, new_allowable_fns)) if '{' and '}' in line else f.write(line)
    if not quiet:
      print(f'|>{(round(36*((file_number+1)/total_files))*"≈" + ">"):<37s}|', end='\r')

def snake(args, quiet=False):
  not_evaluated, permutations, key_file, new_constants, new_allowable_fns = {}, {}, None, copy(constants), copy(allowable_fns)
  # Read external functions from specified file (if callable)
  if args.library:
    path.append(dirname(realpath(args.library)))
    new_allowable_fns.update({k : f for k, f in vars(import_module(psplit(args.library)[-1].rsplit('.')[0])).items() if callable(f) and not f.__name__.startswith("__")})
  # Read the lines in the file
  with open(args.input, 'r') as f:
    input_lines = f.readlines()
  # Read in all constants and lists of perumtation values
  [process_line(line, permutations, not_evaluated, new_constants, new_allowable_fns) for line in input_lines if '{' and '}' in line]
  # Attempt to evaluate any expressions in keys that could not be evaluated on first pass (save expressions that are still unresolved, i.e., those depending on key values)
  resolved = []
  for k, v in not_evaluated.items():
    try:
      result = eval(v, new_allowable_fns, new_constants)
      if result is not None:
        permutations[k] = result
        resolved.append(k)
      else:
        continue
    except:
      pass
  not_evaluated = {k : v for k, v in not_evaluated.items() if k not in resolved}
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
  for fnum, t in enumerate(product(*list(permutations.values()))):
    # Copy existing constants
    iter_constants = copy(new_constants)
    # Get the indicies of the permutation
    idxs = [str(where(asarray(permutations[k]) == t[i])[0][0]) for i, k in enumerate(permutations)]
    # Add current permutation to variable dict
    for k, v in zip(permutations, t):
      iter_constants[k] = v
    # (Re)evaluate all expressions that depend on key values and update evaluator with new variables
    while any(v is None for v in iter_constants.values()):
      num_none = sum(v is None for v in iter_constants.values())
      for k, v in not_evaluated.items():
        iter_constants[k] = evaluate(k, v, iter_constants, not_evaluated, new_allowable_fns)
      if sum(v is None for v in iter_constants.values()) == num_none:
        if not quiet:
          print('The following variables and expressions evaluate to None:')
          for name, none_expr in [(k, v) for k, v in not_evaluated.items() if evaluate(k, v, iter_constants, not_evaluated, new_allowable_fns) is None]:
            print(f' - {name} : {none_expr}')
        raise RecursionError('Expressions are not being resolved! Ensure that user-defined functions cannot return None.')
    # Construct file name based on input arguments
    file_name = get_file_name(args, idxs, t)
    # Write current permutation to file
    write_file(args, file_name, input_lines, fnum, nfiles, iter_constants, not_evaluated, new_allowable_fns, quiet)
    if args.keyfile:
      key_file.write(f'{file_name} {",".join(str(tv) for tv in t)}{linesep}')
  if args.keyfile:
    key_file.close()
  if not quiet:
    print(f'{linesep*2}All {nfiles} files successfully created.')

if __name__ == '__main__':
  snake(get_arguments())
