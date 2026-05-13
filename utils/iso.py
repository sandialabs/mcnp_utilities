#!/usr/bin/env python3

from math import isclose
from argparse import ArgumentParser, RawDescriptionHelpFormatter
# Local modules
from mcnp_utilities.lib.materials import Material, mix_materials, get_compendium_material


def lookup_material(lib, identifier):
  if identifier.isnumeric():
    return get_compendium_material(lib, index=int(identifier))
  else:
    return get_compendium_material(lib, name=identifier)

def determine_components(s):
  level, components, start = 0, [], 0
  for current, c in enumerate(s):
    if c == '(' or c == '{':
      level += 1
    elif c == ')' or c == '}':
      level -= 1

    if (c == ',' or current == len(s) - 1) and level == 0:
      components.append(s[start:current+1].rstrip(',').lstrip(','))
      start = current
  return components

def process_input(inp_str, lvl):
  comps = determine_components(inp_str[1:-1] if (inp_str.startswith('(') and inp_str.endswith(')')) else inp_str)
  if len(comps) > 1:
    fracs = [float(comp.rsplit(':', maxsplit=1)[1]) for comp in comps]
    lvl += 1
    if not (all([f > 0 for f in fracs]) or all([f < 0 for f in fracs])):
      raise ValueError(f'Not all component fractions at level {lvl} have the same sign!')
    if not isclose(abs(sum(fracs)), 1.):
      frac_type = 'Atom' if sum(fracs) > 0 else 'Weight'
      print(f'Note: {frac_type} fractions for component {inp_str} at depth level {lvl} do not sum to unity. Sum: {abs(sum(fracs)):g}')
    return mix_materials([process_input(comp.rsplit(':', maxsplit=1)[0], lvl) for comp in comps], [abs(f) for f in fracs], 'atom' if all([f > 0 for f in fracs]) else 'weight')
  else:
    lvl -= 1
    if ':' in comps[0]:
      comps[0] = comps[0].split(':')[0].lstrip('(').rstrip(')')
    return lookup_material(comps[0].split('=')[0], comps[0].split('=')[1]) if '=' in comps[0] else Material('atom', nuclides={comps[0] : 1})

def parse_arguments():
  parser = ArgumentParser(description=r'''

    _|_|_|    _|_|_|    _|_|
      _|    _|        _|    _|
      _|      _|_|    _|    _|
      _|          _|  _|    _|
    _|_|_|  _|_|_|      _|_|

    (Isotopic Structure Organizer)
    A tool to process and print materials for MCNP''', formatter_class=RawDescriptionHelpFormatter)
  split_combine = parser.add_mutually_exclusive_group()
  parser.add_argument(
    'zaid_pairs',
    type=str,
    help='Parentheses-nested material specifications, with atom or weight fractions specified after the material. Each material consists of a comma-delimited list of ZAID/fraction pairs, with each pair formatted as "<isotope>:<fraction>". "isotope" can be any one of ZAID (ZZAAA), element symbol (Sy[-AAA]), or element name (Name[-AAA]). If AAA is omitted from an entry, or is equal to 0, the entry is considered to be composed of the naturally-occuring isotopes for that element. Atom fractions are positive and weight fractions are negative. All fraction types must match for a given mixture (but do not have to be the same for all mixtures). A comma-delimited list (inside curly braces) of isotopes to exclude when splitting by natural abundance can be included after the isotope name to facilitate the creation of a material enriched with a certain isotope. E.g. U-235:-0.2,U{235}:-0.8. Tip: You may need to enclose this argument in quotes if it contains curly braces.',
    metavar='<ZAID-1[{exZAID-1+...}]:fraction-1>[,<ZAID-2:fraction-2>...]'
  )
  split_combine.add_argument(
    '-s', dest='split_elements',
    action='store_true',
    help='Split elements into isotopic constituents based on natural abundance'
  )
  parser.add_argument(
    '-ec', dest='elemental_carbon',
    action='store_true',
    help='Do not split elemental Carbon into isotopic constituents (for compatibility with pre-ENDF/B VIII.0 libraries)'
  )
  parser.add_argument(
    '-m', dest='print_mat_card',
    nargs=2,
    help='Print formatted material card for MCNP rather than a summary table. Provide material number and fraction type ("atom" or "weight") as whitespace-delimited entries.',
    metavar='<material #> <atom|weight>'
  )
  parser.add_argument(
    '-f', dest='mat_file',
    nargs=3,
    help='Destination file to which to write formatted material card, material number, and fraction type as whitespace-delimited entries.',
    metavar=('<file path>', '<material #>', '<atom|weight>'),
  )
  parser.add_argument(
    '-t', dest='table',
    nargs=2,
    help='Print material specification (in desired column order as a comma-delimited list) as a LaTeX table to a file (Z=ZAID,N=Name,A=Atom Fraction,W=Weight Fraction). Format specifiers can optionally be appended to A or W entries after a colon (e.g. "W:.4f"). The default format specifier is ".6E".',
    metavar=('<file path>', '<format specifiers>')
  )
  parser.add_argument(
    '-d', dest='density',
    type=float,
    help='Atom (positive) [a/(b-cm)] or mass (negative) [g/cc] density of the material. If provided, a breakdown of component atom densities will be provided in the summary table.',
    metavar='<density>'
  )
  parser.add_argument(
    '-n', dest='name',
    type=str,
    help='Name of the material (optional)',
    metavar='<name>'
  )
  split_combine.add_argument(
    '-c', dest='combine',
    action='store_true',
    help='Combine isotopic constituents into elemental components.'
  )
  return parser.parse_args()

if __name__ == '__main__':
  args = parse_arguments()
  mat_level = 0
  mat = process_input(args.zaid_pairs, mat_level)

  if args.split_elements:
    mat.split_elements(split_carbon=not args.elemental_carbon)
  elif args.combine:
    mat.combine_isotopes()

  if args.name:
    mat.name = args.name

  if args.density:
    mat.density = args.density

  if args.print_mat_card or args.mat_file:
    if args.print_mat_card:
      mat.number, mat_type = args.print_mat_card
      print(mat.to_MCNP_material_card(mat_type))
    if args.mat_file:
      fpath, mat.number, mat_type = args.mat_file
      with open(fpath, 'w') as f:
        f.write(mat.to_MCNP_material_card(mat_type))
  else:
    mat.print_table()

  if args.table:
    mat.write_latex_table(args.table[0], [c for c in args.table[1].split(',')])

  print(f'Mixture molar mass (per atom): {mat.get_molar_mass():g} g/mol')
