#!/usr/bin/env python3

# Author: Zackary Dodson


from argparse import ArgumentParser, RawDescriptionHelpFormatter
from bisect import bisect_left
from collections import OrderedDict
from copy import copy
from itertools import product
from multiprocessing import Pool, cpu_count, freeze_support
from math import pi, ceil
from numpy import zeros, empty, asarray, cos, sin, size, ascontiguousarray, shape
from pyevtk.hl import unstructuredGridToVTK
from pyevtk.vtk import VtkVoxel
from h5py import File as h5file


def write_pvtu_file(filelist, name, setname, Es):
  namen = f'{name}.pvtu'
  Es = [f'_{e}' for e in Es] if len(Es) > 1 else ['']
  with open(namen, 'w') as f:
    f.write('<?xml version="1.0"?>\n')
    f.write('<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
    f.write('  <PUnstructuredGrid GhostLevel="0">\n')
    f.write('    <PPoints>\n')
    f.write('      <PDataArray type="Float64" NumberOfComponents="3" format="appended"/>\n')
    f.write('    </PPoints>\n')
    f.write('    <PCells>\n')
    f.write('      <DataArray type="Int32" Name="connectivity" format="appended"/>\n')
    f.write('      <DataArray type="Int32" Name="offsets" format="appended"/>\n')
    f.write('      <DataArray type="Int32" Name="types" format="appended"/>\n')
    f.write('    </PCells>\n')
    f.write('    <PCellData>\n')
    for e in Es:
      f.write(f'      <DataArray type="Float64" Name="{setname}{e}" format="appended"/>\n')
      f.write(f'      <DataArray type="Float64" Name="{setname}_rel_err{e}" format="appended"/>\n')
    f.write('    </PCellData>\n')
    [f.write(f'    <Piece Source="{f}.vtu"/>\n') for f in filelist]
    f.write('  </PUnstructuredGrid>\n')
    f.write('</VTKFile>')

def find_bounds(val, vals):
  idx = bisect_left(vals, val)
  return [vals[idx - 1], vals[idx]]

def procData(data):
  pts, data, ncell, xs, ys, zs, name, dataname, is_cyl, nEmesh, Emesh = data[0:11]
  cellpoints, cell_list, n_verts = empty([8, 3]), [], ncell * 8
  vs = [xs, ys, zs]
  for n in range(ncell):
    cellpoints[:, :] = list(product(*[find_bounds(pts[n, i], vs[i]) for i in range(3)]))
    c = cell(copy(cellpoints), data[n, 0, :], data[n, 1, :])
    cell_list.append(c)
  if is_cyl:
    [c.convert_to_rec() for c in cell_list]
  x, y, z, dat, err = empty((n_verts)), empty((n_verts)), empty((n_verts)), empty((ncell, nEmesh)), empty((ncell, nEmesh))
  for j, cel in enumerate(cell_list):
    dat[j, :], err[j, :] = cel.values, cel.rel_errors
    for i in range(8):
      idx = j * 8 + i
      x[idx], y[idx], z[idx] = cel.points[i, :]
  cellData = OrderedDict()
  for e in range(nEmesh):
    E = f'_{Emesh[e]}' if nEmesh > 1 else ''
    cellData[f'FMESH{dataname}{E}'] = ascontiguousarray(dat[:, e])
    cellData[f'FMESH{dataname}_rel_err{E}'] = ascontiguousarray(err[:, e])
  unstructuredGridToVTK(
    name,
    x,
    y,
    z,
    asarray([v for v in range(n_verts)]),
    asarray([v * 8 for v in range(1, ncell + 1)]),
    asarray([VtkVoxel.tid for v in range(ncell)]),
    cellData
  )

def proc_vtu_data(Dset, name, nfiles, tasks):
  inplist, namelist = [], []
  inc = upp = ceil(Dset.ncells / nfiles)
  dsetname, Emesh = Dset.name, Dset.Emesh
  for i in range(nfiles):
    namen = f'{name}_FMESH{Dset.name}_{i}' if nfiles > 1 else f'{name}_FMESH{Dset.name}'
    namelist.append(namen)
    inplist.append(
      [
        Dset.mid_cell_pts[(i * inc):upp, :],
        Dset.cellData[(i * inc):upp, :, :],
        upp - (i * inc),
        Dset.Xbounds,
        Dset.Ybounds,
        Dset.Zbounds,
        namen,
        Dset.name,
        Dset.cyl,
        Dset.nEmesh,
        Dset.Emesh
      ]
    )
    upp += inc
    upp = Dset.ncells if upp >= Dset.ncells else upp
  del Dset
  with Pool(processes=tasks) as pool:
    pool.map(procData, inplist)
  if nfiles > 1:
    write_pvtu_file(namelist, f'{name}_FMESH{dsetname}', f'FMESH{dsetname}', Emesh)

class cell:
  def __init__(self, p, z, e):
    self.points, self.values, self.rel_errors = p, z, e

  def convert_to_rec(self):
    '''
    Note: Theta is in units of revolutions.
    '''
    temp = zeros([8, 3])
    temp[:, 0] = self.convert_r_to_x(self.points[:, 0], self.points[:, 2] * 2 * pi)
    temp[:, 1] = self.convert_z_to_y(self.points[:, 0], self.points[:, 2] * 2 * pi)
    temp[:, 2] = self.points[:, 1]
    self.points = temp

  def convert_r_to_x(self, r_value, theta_value):
    '''
    Returns the x value or vector corresponding to the input radius and
    theta values (or vectors).

    Note: Input theta should be in units of radians.
    '''
    return r_value * cos(theta_value)

  def convert_z_to_y(self, r_value, theta_value):
    '''
    Returns the y value or vector corresponding to the input radius and
    theta values (or vectors).

    Note: Input theta should be in units of radians.
    '''
    return r_value * sin(theta_value)

class Dataset:
  def __init__(self, xpts, ypts, zpts, xm, ym, zm, data, name, cylindical=False):
    self.Xbounds, self.Ybounds, self.Zbounds = xpts, ypts, zpts
    self.x_mid, self.y_mid, self.z_mid = xm, ym, zm
    self.data, self.name, self.cyl = data, name, cylindical
    self.parsed = False
    self.parse_data()

  def parse_data(self):
    self.Emesh = list(self.data.keys())
    self.nEmesh, self.ncells = len(self.Emesh), size(self.data[self.Emesh[0]], 0)
    self.cellData, self.mid_cell_pts = empty((self.ncells, 2, self.nEmesh)), empty((self.ncells, 3))
    for i, p in enumerate(product(self.x_mid, self.y_mid, self.z_mid)):
      self.mid_cell_pts[i,:] = p
    for i, e in enumerate(self.Emesh):
      self.cellData[:, :, i] = self.data[e][:, :]
    del self.data

def read_h5_datasets(h5g):
  Dsets = []
  if 'results' in h5g:
    for tally_p in h5g['results/mesh_tally']:
      tnum = tally_p.split('_')[-1]
      mesh_path = f'results/mesh_tally/mesh_tally_{tnum}'
      h5gd = h5g[mesh_path]
      e_bounds, A, CYLflag = asarray(h5gd['grid_energy'][1:]), {}, 'grid_r' in h5gd
      l = ['r', 'z', 't'] if CYLflag else ['x', 'y', 'z']
      X, Y, Z = [asarray([(h5gd[f'grid_{d}'][i]+h5gd[f'grid_{d}'][i+1])/2 for i in range(len(h5gd[f'grid_{d}'])-1)]) for d in l]
      X_b, Y_b, Z_b = (asarray(h5gd[f'grid_{d}']) for d in l)
      data, err = asarray(h5gd['mean'][:, 0, :, :, :]).T, asarray(h5gd['relative_standard_error'][:, 0, :, :, :]).T
      sd1, sd2, sd3, nEmesh = shape(data)
      ncells = sd1*sd2*sd3
      for e in range(nEmesh):
        e_ub = str(e_bounds[e])
        A[e_ub] = empty((ncells, 2))
        for idx, (i, j, k) in enumerate(product(range(sd1), range(sd2), range(sd3))):
          A[e_ub][idx,:] = data[i,j,k,e], err[i,j,k,e]
      Dsets.append(Dataset(X_b, Y_b, Z_b, X, Y, Z, A, tnum, CYLflag))
  else:
    for tnum in h5g:
      h5gd = h5g[tnum]
      e_bounds, A, CYLflag = asarray(h5gd['E_bounds']), {}, 'T_bounds' in h5gd
      if CYLflag:
        X, Y, Z = (asarray(h5gd[d]) for d in ['R', 'Z', 'T'])
        X_b, Y_b, Z_b = (asarray(h5gd[d]) for d in ['R_bounds', 'Z_bounds', 'T_bounds'])
      else:
        X, Y, Z = (asarray(h5gd[d]) for d in ['X', 'Y', 'Z'])
        X_b, Y_b, Z_b = (asarray(h5gd[d]) for d in ['X_bounds', 'Y_bounds', 'Z_bounds'])
      data, err = asarray(h5gd['data']), asarray(h5gd['rel_err'])
      sd1, sd2, sd3, nEmesh = shape(data)
      ncells = sd1*sd2*sd3
      for e in range(nEmesh):
        e_ub = str(e_bounds[e])
        A[e_ub] = empty((ncells, 2))
        for idx, (i, j, k) in enumerate(product(range(sd1), range(sd2), range(sd3))):
          A[e_ub][idx,:] = data[i,j,k,e], err[i,j,k,e]
      if ('total_data' in h5gd) and ('total_rel_err' in h5gd):
        total_data, total_err = asarray(h5gd['total_data']), asarray(h5gd['total_rel_err'])
        A['Total'] = empty((ncells, 2))
        for idx, (i, j, k) in enumerate(product(range(sd1), range(sd2), range(sd3))):
          A['Total'][idx,:] = total_data[i,j,k], total_err[i,j,k]
      Dsets.append(Dataset(X_b, Y_b, Z_b, X, Y, Z, A, tnum, CYLflag))
  return Dsets

def get_h5_datasets(fname):
  with h5file(fname, 'r') as f:
    return read_h5_datasets(f)

def run_paramesh(name, tasks, nfiles):
  [proc_vtu_data(D, name, nfiles, tasks) for D in get_h5_datasets(name)]

def get_arguments():
  parser = ArgumentParser(
    description='''
 ██████╗  █████╗ ██████╗  █████╗ ███╗   ███╗███████╗███████╗██╗  ██╗
 ██╔══██╗██╔══██╗██╔══██╗██╔══██╗████╗ ████║██╔════╝██╔════╝██║  ██║
 ██████╔╝███████║██████╔╝███████║██╔████╔██║█████╗  ███████╗███████║
 ██╔═══╝ ██╔══██║██╔══██╗██╔══██║██║╚██╔╝██║██╔══╝  ╚════██║██╔══██║
 ██║     ██║  ██║██║  ██║██║  ██║██║ ╚═╝ ██║███████╗███████║██║  ██║
 ╚═╝     ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝╚══════╝╚══════╝╚═╝  ╚═╝
    A tool to convert MCNP mesh tally files to vtu files.''',
    formatter_class=RawDescriptionHelpFormatter
  )
  parser.add_argument(
    '-j', '--processes',
    help='Number of threads to use to process files. (default: %(default)s)',
    default=1,
    type=int,
    metavar='<# processes>'
  )
  parser.add_argument(
    '-n', '--numfiles',
    help='Number of files to split vtu files into. (default: %(default)s)',
    default=1,
    type=int,
    metavar='<# output files>'
  )
  parser.add_argument(
    'files',
    help='Files to process (whitespace delimited list) (only supports hdf5).',
    nargs='+',
    metavar='<file>'
  )
  return parser.parse_args()

def main(args):
  if (args.processes > args.numfiles):
    print(f'Only {args.numfiles} process(es) will be utilized.\n')
    args.processes = args.numfiles
  if args.processes > 1:
    totprocesses = cpu_count()
    if args.processes > totprocesses:
      raise RuntimeError(f'Error, only {totprocesses} processes are available. Terminating.')
  num_files = len(args.files)
  for i, f in enumerate(args.files, start=1):
    print(f'Processing file "{f}" ({i}/{num_files})...')
    run_paramesh(f, args.processes, args.numfiles)

if __name__ == '__main__':
  freeze_support()
  main(get_arguments())
