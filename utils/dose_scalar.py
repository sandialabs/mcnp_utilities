#!/usr/bin/env python3

from pickle import load
from argparse import ArgumentParser, Action
from os.path import join as pjoin, dirname
# Local modules
from mcnp_utilities.utils.iso import process_input


mass_db = pjoin(dirname(__file__), '..', 'data', 'iso_mass.pkl')     # Mass database file (taken from xsdir file)

N_A = 5.9703109E+23                         # Avogadro's Number (neutron relative) [a/mol]
s_in_hr = 3600                              # [s/hr]
MeVg_to_rad = 1. / (100. * 624150.64799632) # [rad/(MeV/g)]
cm2_in_b = 1.0E-24                          # [cm^2/b]

nu_bar = 2.44                               # [n/f]
Q_bar = 192.4                               # [MeV/f]
MeV_to_MJ = 1.602176487E-19                 # [MJ/MeV]
C = nu_bar / (Q_bar * MeV_to_MJ)            # [n/MJ]

def load_mass_database(mass_datafile):
  with open(mass_datafile, 'rb') as f:
    return load(f)

class CommaSplitAction(Action):
  def __call__(self, parser, namespace, values, option_string=None):
    setattr(namespace, self.dest, values.split(","))

def get_arguments():
  parser = ArgumentParser(description='A tool to compute the FM card scalar for an MCNP material to obtain dose [rad/(MW-hr)] in that material.')
  parser.add_argument(
    'zaid_pairs',
    type=str,
    help='Parentheses-nested material specifications, with atom or weight fractions specified after the material. Each material consists of a comma-delimited list of ZAID/fraction pairs, with each pair formatted as "<isotope>:<fraction>". "isotope" can be any one of ZAID (ZZAAA), element symbol (Sy[-AAA]), or element name (Name[-AAA]). If AAA is omitted from an entry, or is equal to 0, the entry is considered to be composed of the naturally-occuring isotopes for that element. Atom fractions are positive and weight fractions are negative. All fraction types must match for a given mixture (but do not have to be the same for all mixtures). A comma-delimited list (inside curly braces) of isotopes to exclude when splitting by natural abundance can be included after the isotope name to facilitate the creation of a material enriched with a certain isotope. E.g. U-235:-0.2,U{235}:-0.8. Tip: You may need to enclose this argument in quotes if it contains curly braces.',
    metavar='<ZAID-1[{exZAID-1+...}]:fraction-1>[,<ZAID-2:fraction-2>...]'
  )
  return parser.parse_args()


if __name__ == '__main__':
  # Process arguments
  args = get_arguments()
  # Load isotope molar masses (relative to neutron mass)
  db = load_mass_database(mass_db)
  # Create material from input
  mat = process_input(args.zaid_pairs)
  mat.split_elements()
  # Compute the molar mass of the mixture
  M_mix = sum([af * db[iso] for iso, af in mat.nuclide_atom_fractions.items()])
  # Print the FM muliplier for the material to obtain dose in rad/(MW-hr)
  print(f'{N_A * s_in_hr * C * MeVg_to_rad * cm2_in_b / M_mix:.6E} [rad[material]/(MW*hr)]')
  print(f'{N_A * C * MeVg_to_rad * cm2_in_b / M_mix:.6E} [rad[material]/MJ]')
