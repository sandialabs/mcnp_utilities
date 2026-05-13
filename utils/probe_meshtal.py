#!/usr/bin/env python3

from h5py import File as h5file
from uncertainties import ufloat, UFloat
from argparse import ArgumentParser
from numpy import asarray, pow, sqrt


def quad_sum(arr, axs):
  return sqrt(pow(arr, 2).sum(axis=axs))

def get_arguments():
  parser = ArgumentParser(description='A tool to probe meshtal data')
  parser.add_argument(
    'file',
    type=str
  )
  parser.add_argument(
    'location',
    help='Location to probe',
    type=str
  )
  parser.add_argument(
    '-s', dest='scalar',
    type=float,
    default=1.0
  )
  return parser.parse_args()

def meshtal_probe(meshtal_file, tally, iloc, jloc, kloc, eloc=None) -> UFloat:
  tally = str(tally)
  with h5file(meshtal_file) as f:
    if 'X' and 'Y' and 'Z' in f[tally]:
      x_lb, x_ub = asarray(f[tally]['X_bounds'])[0], asarray(f[tally]['X_bounds'])[-1]
      y_lb, y_ub = asarray(f[tally]['Y_bounds'])[0], asarray(f[tally]['Y_bounds'])[-1]
      z_lb, z_ub = asarray(f[tally]['Z_bounds'])[0], asarray(f[tally]['Z_bounds'])[-1]

      if not (x_lb <= iloc <= x_ub):
        print(f'Warning: specified x location of {iloc} cm is outside of x-bounds ({x_lb} - {x_ub}) cm!')
      if not (y_lb <= jloc <= y_ub):
        print(f'Warning: specified y location of {jloc} cm is outside of y-bounds ({y_lb} - {y_ub}) cm!')
      if not (z_lb <= kloc <= z_ub):
        print(f'Warning: specified z location of {kloc} cm is outside of z-bounds ({z_lb} - {z_ub}) cm!')

      i, j, k = (abs(asarray(f[tally]['X']) - iloc)).argmin(), (abs(asarray(f[tally]['Y']) - jloc)).argmin(), (abs(asarray(f[tally]['Z']) - kloc)).argmin()
    elif 'R' and 'Z' and 'T' in f[tally]:
      r_lb, r_ub = asarray(f[tally]['R_bounds'])[0], asarray(f[tally]['R_bounds'])[-1]
      z_lb, z_ub = asarray(f[tally]['Z_bounds'])[0], asarray(f[tally]['Z_bounds'])[-1]
      t_lb, t_ub = asarray(f[tally]['T_bounds'])[0], asarray(f[tally]['T_bounds'])[-1]

      if not (r_lb <= iloc <= r_ub):
        print(f'Warning: specified r location of {iloc} cm is outside of r-bounds ({r_lb} - {r_ub}) cm!')
      if not (z_lb <= jloc <= z_ub):
        print(f'Warning: specified z location of {jloc} cm is outside of z-bounds ({z_lb} - {z_ub}) cm!')
      if not (t_lb <= kloc <= t_ub):
        print(f'Warning: specified t location of {kloc} revolution is outside of t-bounds ({t_lb} - {t_ub}) revolutions!')

      i, j, k = (abs(asarray(f[tally]['R']) - iloc)).argmin(), (abs(asarray(f[tally]['Z']) - jloc)).argmin(), (abs(asarray(f[tally]['T']) - kloc)).argmin()

    if eloc is not None:
      for e_idx, e_ub in enumerate(f[tally]['E_bounds']):
        if eloc < e_ub:
          e = e_idx
      data_value, rel_err = asarray(f[tally]['data'])[i, j, k, e], asarray(f[tally]['rel_err'])[i, j, k, e]
    else:
      data_value, rel_err = asarray(f[tally]['data']).sum(axis=3)[i, j, k], quad_sum(asarray(f[tally]['rel_err']), 3)[i, j, k]

  return ufloat(data_value, data_value * rel_err)

if __name__ == '__main__':
  args = get_arguments()
  locs = args.location.split(',')
  if len(locs) == 4:
    iloc, jloc, kloc, eloc = [float(v) for v in locs]
  elif len(locs) == 3:
    iloc, jloc, kloc = [float(v) for v in locs]
  fyle, tally_number = args.file.split(':')

  if len(locs) == 4:
    print(f'{meshtal_probe(fyle, tally_number, iloc, jloc, kloc, eloc=eloc) * args.scalar:.1uS}')
  elif len(locs) == 3:
    print(f'{meshtal_probe(fyle, tally_number, iloc, jloc, kloc) * args.scalar:.1uS}')
