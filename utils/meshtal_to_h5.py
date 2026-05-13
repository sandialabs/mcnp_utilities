#!/usr/bin/env python3

# Author: Zackary Dodson

from itertools import product
from numpy import empty, asarray
from h5py import File as h5file
from argparse import ArgumentParser


def get_tallies(fname):
  nums = []
  with open(fname, 'r') as f:
    for line in f:
      line = line.split()
      if len(line) == 4:
        if line[:3] == ["Mesh", "Tally", "Number"]:
          nums.append(line[3])
  return nums

def get_tally_lines(fname, tally_num):
  tally_lines, in_tally = [], False
  with open(fname, 'r') as f:
    for line in f:
      line = line.split()
      if line:
        if line[:3] == ["Mesh", "Tally", "Number"]:
          in_tally = True if line[3] == tally_num else False
        if in_tally:
          tally_lines.append(line)
  return tally_lines

def read_data(line_list,dim_tuple):
  data, err, return_dict = empty(dim_tuple), empty(dim_tuple), {}
  for idx, line in enumerate(line_list):
    if all([x in line for x in ['Result', 'Error']]):
      res_idx = line.index('Result')
      for offset, (e, i, j, k) in enumerate(product(*(range(dim_tuple[d]) for d in [3, 0, 1, 2]))):
        lidx = idx + offset + 1
        data[i,j,k,e], err[i,j,k,e] = float(line_list[lidx][res_idx]), float(line_list[lidx][res_idx+1])
      if dim_tuple[3] > 1:
        total_data, total_err = empty(dim_tuple[:-1]), empty(dim_tuple[:-1])
        for offset, (i, j, k) in enumerate(product(*(range(d) for d in dim_tuple[:3]))):
          lidx2 = lidx + offset + 1
          total_data[i,j,k], total_err[i,j,k] = float(line_list[lidx2][res_idx]), float(line_list[lidx2][res_idx+1])
  return_dict['data'], return_dict['rel_err'], return_dict['comment'] = data, err, " ".join(line_list[1])
  if dim_tuple[3] > 1:
    return_dict['total_data'], return_dict['total_rel_err'] = total_data, total_err
  return return_dict

def generate_xyz_dataset(line_list):
  d4len = 1
  for line in line_list:
    if line[-1] != 'Error':
      if line[0] == 'X':
        d1len, d1bounds = len(line) - 3, asarray(line[2:], dtype=float)
      elif line[0] == 'Y':
        d2len, d2bounds = len(line) - 3, asarray(line[2:], dtype=float)
      elif line[0] == 'Z':
        d3len, d3bounds = len(line) - 3, asarray(line[2:], dtype=float)
      elif line[:2] == ['Energy', 'bin']:
        d4len, d4bounds = len(line) - 4, asarray(line[4:], dtype=float)
    else:
      break
  dat = read_data(line_list,(d1len, d2len, d3len, d4len))
  dat['X_bounds'], dat['Y_bounds'], dat['Z_bounds'], dat['E_bounds'] = d1bounds, d2bounds, d3bounds, d4bounds
  for d, a in zip(['X', 'Y', 'Z'], [d1bounds, d2bounds, d3bounds]):
    dat[d] = asarray([(a[i]+a[i+1])/2 for i in range(len(a)-1)])
  dat['dims'] = ['X(cm)','Y(cm)','Z(cm)','Energy(MeV)']
  return dat

def generate_cyl_dataset(line_list):
  d4len = 1
  for line in line_list:
    if line[-1] != 'Error':
      if line[:2] == ['R', 'direction:']:
        d1len, d1bounds = len(line) - 3, asarray(line[2:], dtype=float)
      elif line[0] == 'Z':
        d2len, d2bounds = len(line) - 3, asarray(line[2:], dtype=float)
      elif line[0] == 'Theta':
        d3len, d3bounds = len(line) - 4, asarray(line[3:], dtype=float)
      elif line[:2] == ['Energy', 'bin']:
        d4len, d4bounds = len(line) - 4, asarray(line[4:], dtype=float)
    else:
      break
  dat = read_data(line_list,(d1len, d2len, d3len, d4len))
  dat['R_bounds'], dat['Z_bounds'], dat['T_bounds'], dat['E_bounds'] = d1bounds, d2bounds, d3bounds, d4bounds
  for d, a in zip(['R', 'Z', 'T'], [d1bounds, d2bounds, d3bounds]):
    dat[d] = asarray([(a[i] + a[i+1])/2 for i in range(len(a)-1)])
  dat['dims'] = ['R(cm)','Z(cm)','Theta(rev)','Energy(MeV)']
  return dat

def parse_tally(line_list):
  dirs = []
  for line in line_list:
    if line[1].startswith('direction'):
      dirs.append(line[0])
  return generate_cyl_dataset(line_list) if 'Theta' in dirs else generate_xyz_dataset(line_list)

def get_dataset(fname, tally_num):
  return parse_tally(get_tally_lines(fname, tally_num))

def write_file(h5_name,meshtal_name,procs=1):
  if procs > 1:
    from multiprocessing import Pool
    tals = get_tallies(meshtal_name)
    with Pool(processes=procs) as pool, h5file(h5_name, 'w') as f:
      dsets = {t : d for t, d in zip(tals, pool.starmap(get_dataset, [(meshtal_name,t) for t in tals]))}
      for t, dset in dsets.items():
        tal_dset = f.create_group(t)
        for k, v in dset.items():
          tal_dset[k] = v
  else:
    with h5file(h5_name, 'w') as f:
      for t in get_tallies(meshtal_name):
        tal_dset = f.create_group(t)
        res = get_dataset(meshtal_name,t)
        for k, v in res.items():
          tal_dset[k] = v

def get_arguments():
  parser = ArgumentParser(description="A tool to convert an MCNP meshtal file to an hdf5 file.",)
  parser.add_argument(
    '-o', '--output',
    help='Name of the output hdf5 file. Default: <input>.h5)',
    metavar='<output file>'
  )
  parser.add_argument(
    'file',
    help='MCNP meshtal file to process.',
    metavar='<input>'
  )
  parser.add_argument(
    '-j', '--processes',
    type=int,
    help='Number of processes to use (speeds up conversion of files with multiple tallies). Default=%(default)s',
    metavar='# processes',
    default=1
  )
  return parser,parser.parse_args()

if __name__ == '__main__':
  aparser, args = get_arguments()
  args.output = f'{args.file}.h5' if not args.output else args.output
  write_file(args.output,args.file,procs=args.processes)
