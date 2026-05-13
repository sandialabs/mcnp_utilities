#! /usr/bin/env python3

from argparse import ArgumentParser, Action
from numpy import asarray, empty
from math import log
from pickle import dump
from uncertainties import ufloat, UFloat
# Local modules
from mcnp_utilities.utils.tallies import get_all_tallies
from mcnp_utilities.utils.get_rx import get_keff


class SpectrumTally:
  def __init__(self, energy_bin_edges, tally_data, sdev, keff, Emin=1.E-50):
    self.energy_bin_edges = asarray(energy_bin_edges).astype(float)
    data = asarray(tally_data).astype(float)
    sdev = asarray(sdev).astype(float)
    self.tally_data = empty((len(data)), dtype=UFloat.dtype)
    for i, (v, e) in enumerate(zip(data, sdev)):
      self.tally_data[i] = ufloat(v, e)
    self.keff = keff
    self.bin_number = len(self.tally_data)
    self.bin_widths, self.u, self.bin_midpoints = [], [], []
    for i in range(len(self.energy_bin_edges)):
      if i > 0:
        self.bin_widths.append(self.energy_bin_edges[i] - self.energy_bin_edges[i-1])
        self.u.append(log(self.energy_bin_edges[i] / self.energy_bin_edges[i-1]))
        self.bin_midpoints.append((self.energy_bin_edges[i] + self.energy_bin_edges[i-1]) / 2)
      else: # fill the lowest bin lethargy approximately i.e. Emin = 1e-50
        self.bin_widths.append(self.energy_bin_edges[i] - Emin)
        self.u.append(log(self.energy_bin_edges[i] / Emin))
        self.bin_midpoints.append((self.energy_bin_edges[i] + Emin) / 2)
    self.bin_widths = asarray(self.bin_widths).astype(float)
    self.u = asarray(self.u).astype(float)
    self.tally_data /= (self.u * self.keff) # Divide by lethargy and keff
    self.bin_midpoints = asarray(self.bin_midpoints).astype(float)

class ColonSplitAction(Action):
  def __call__(self, parser, namespace, values, option_string=None):
    setattr(namespace, self.dest, values.split(":"))

def get_arguments():
  parser = ArgumentParser(description='A tool to extract energy spectra from MCNP output files')
  parser.add_argument(
    'tally_path',
    type=str,
    action=ColonSplitAction,
    help='Path to tally delimited by colons',
    metavar='<file>:<tally #>:<cell #>:<material #>'
  )
  parser.add_argument(
    'output_file',
    type=str,
    help='Name of the file to which to write tally information',
    metavar='<output file name>'
  )
  parser.add_argument(
    '-m', dest='manual',
    action='store_true',
    help='Flag to set input file as manual'
  )
  return parser.parse_args()

if __name__ == '__main__':
  args = get_arguments()
  with open(args.tally_path[0], 'r', errors='ignore') as fout, open(args.output_file, 'wb') as fspec:
    if not args.manual:
      tally = get_all_tallies(fout.readlines())[args.tally_path[1]][args.tally_path[2]][args.tally_path[3]]
      dump(
        SpectrumTally(
          [k for k in tally if k != 'total'],
          [v.nominal_value for k, v in tally.items() if k != 'total'],
          [v.std_dev for k, v in tally.items() if k != 'total'],
          get_keff(args.tally_path[0]).nominal_value if get_keff(args.tally_path[0]) is not None else 1.0
        ),
        fspec
      )
    else:
      energies, data, stddevs = [], [], []
      for line in fout:
        e, d, s = line.split()
        energies.append(e)
        data.append(d)
        stddevs.append(s)
      dump(SpectrumTally(energies, data, stddevs, 1.0), fspec)
