#! /usr/bin/env python3

from numpy import empty, asarray, size, sqrt
from multiprocessing import Pool, freeze_support
from h5py import File as h5file


class Dataset:
  def __init__(self, data, energy_ub, total=None):
    self.data = data
    self.energy_ub = energy_ub
    if total is not None:
      tmp = empty((size(self.data, 0), size(self.data, 1), size(self.data, 2), size(self.data, 3)+1))
      tmp[:,:,:,:-1] = self.data[:,:,:,:]
      tmp[:,:,:,-1] = total[:,:,:]
      self.data = tmp
      self.energy_ub.append('Total')
    if len(energy_ub) == 1 and energy_ub[0] == 1E+36:
      self.energy_ub = ['Total']

def get_h5_dataset(h5name, tal_num):
  with h5file(h5name,'r') as f:
    if 'total_data' in f[tal_num]:
      return Dataset(asarray(f[tal_num]['data']), list(f[tal_num]['E_bounds']), total=asarray(f[tal_num]['total_data']))
    else:
      return Dataset(asarray(f[tal_num]['data']), list(f[tal_num]['E_bounds']))

def rms(arr):
  return sqrt((arr**2).mean())

def get_zeta(arr):
  return rms((arr / arr.mean()) - 1.0)

def process_energy_group(arr, mididx):
  return get_zeta(arr), get_zeta(arr[:,mididx,:]), get_zeta(arr.sum((0,2))), get_zeta(arr.sum(1)), arr.mean(), arr.max()

def process_tally(dset, divisor, procs=1):
  dset.data /= divisor
  midz = size(dset.data, 1) // 2
  with Pool(processes=procs) as pool:
    return {k : v for k, v in zip([dset.energy_ub[e] for e in range(size(dset.data, 3))], pool.starmap(process_energy_group, [(dset.data[:,:,:,e], midz) for e in range(size(dset.data, 3))]))}

def get_arguments():
  parser = ArgumentParser(description="A tool to compute asymmetry metrics (RMS of deviation from mean) for a cylindrical MCNP FMESH tally.")
  parser.add_argument(
    'file',
    help='Relative path to meshtal file (only supports h5 files).',
    type=str,
    metavar='<file>'
  )
  parser.add_argument(
    'tally_number',
    help='Tally number(s) to process (in a comma-delimited list).',
    type=str,
    metavar='<tally number 1>[,tally number 2...]'
  )
  parser.add_argument(
    '-d','--divide-keff',
    help='Divide result by keff in this MCNP out file.',
    type=str,
    metavar='<MCNP output file>'
  )
  parser.add_argument(
    '-j', '--processes',
    help='Number of threads to use to process data. Default: %(default)s',
    type=int,
    default=1,
    metavar='<# processes>'
  )
  return parser.parse_args()

def main():
  args = get_arguments()
  dv = get_keff(args.divide_keff).nominal_value if args.divide_keff else 1.0
  tnums = args.tally_number.split(',')
  for tally_num in tnums:
    print('{0:^110}'.format(tally_num))
    ds = get_h5_dataset(args.file, tally_num)
    res = process_tally(ds, dv, args.processes)
    hdrs = ['EUB (MeV)', 'Asym', 'Asym MP', 'A Asym', 'RT Asym', 'Avg.', 'Max.']
    print('{:^15}|{:^15}|{:^15}|{:^15}|{:^15}|{:^15}|{:^15}'.format(*hdrs))
    print('-'*15 + 5*('+' + '-'*15) + '+' + '-'*14)
    [print('{:^15}|{:^15.5E}|{:^15.5E}|{:^15.5E}|{:^15.5E}|{:^15.5E}|{:^15.5E}'.format(k,*v)) for k, v in res.items()]

if __name__=="__main__":
  from mcnp_utilities.utils.get_rx import get_keff
  from argparse import ArgumentParser

  freeze_support()
  main()
