#!/usr/bin/env python3

from argparse import ArgumentParser, Namespace
from h5py import File as h5file
from matplotlib.pyplot import subplots, ticklabel_format, get_cmap
from matplotlib.colors import LogNorm, ListedColormap, BoundaryNorm
from matplotlib.ticker import MultipleLocator, LogLocator, NullFormatter
from numpy import asarray, abs, nan, pow, sqrt, nanmax, nanmin, pi, vstack, concatenate, arange, ndarray, squeeze
from uncertainties import ufloat
# Local modules
from mcnp_utilities.lib.basic_tools import set_plot_params
from mcnp_utilities.utils.tallies import is_float


units = {
  'mm' : 1e+1,
  'cm' : 1e+0,
  'm'  : 1e-2,
  'km' : 1e-5,
  'ft' : 1./(2.54*12.)
}

class MeshDataset:
  def __init__(self, h5f, tally, e_idx=None):
    t0 = tally.split('+')[0]
    if 'results' in h5f: # Data is in the runtape
      mesh_path = f'results/mesh_tally/mesh_tally_'
      self.is_cyl = all([all(f'{mesh_path}{t}/grid_{c}' in h5f for c in ['r', 'z', 't']) for t in tally.split('+')])
      self.data = sum([asarray(h5f[f'{mesh_path}{t}/mean']) for t in tally.split('+')])
      self.abs_error = sqrt(sum([pow(asarray(h5f[f'{mesh_path}{t}/relative_standard_error']) * asarray(h5f[f'{mesh_path}{t}/mean']), 2) for t in tally.split('+')]))
      e_i = e_idx if e_idx is not None else -1
      self.data, self.abs_error = self.data[e_i, 0, :, :, :].T, self.abs_error[e_i, 0, :, :, :].T
      self.dims, self.dim_bounds, axes = {}, {}, ['r', 'z', 't'] if self.is_cyl else ['x', 'y', 'z']
      for d in axes:
        self.dim_bounds[d.upper()] = asarray(h5f[f'{mesh_path}{t0}/grid_{d}'])
        a = h5f[f'{mesh_path}{t0}/grid_{d.lower()}']
        self.dims[d.upper()] = asarray([(a[i]+a[i+1])/2 for i in range(len(a)-1)])
      if self.is_cyl:
        self.dims['T'] /= 2 * pi
        self.dim_bounds['T'] /= 2 * pi
    else: # Data is in homemade meshtal.h5
      self.is_cyl = all([all(c in h5f[t] for c in ['R', 'Z', 'T']) for t in tally.split('+')])
      self.data = sum([asarray(h5f[t]['data']) for t in tally.split('+')])
      self.abs_error = sqrt(sum([pow(asarray(h5f[t]['rel_err']) * asarray(h5f[t]['data']), 2) for t in tally.split('+')]))
      if e_idx is not None:
        self.data, self.abs_error = self.data[..., e_idx], self.abs_error[..., e_idx]
      else:
        self.data, self.abs_error = self.data.sum(axis=3), sqrt(pow(self.abs_error, 2).sum(axis=3))
      self.dims = {d : asarray(h5f[t0][d]) for d in ['R', 'Z', 'T']} if self.is_cyl else {d : asarray(h5f[t0][d]) for d in ['X', 'Y', 'Z']}
      self.dim_bounds = {d : asarray(h5f[t0][f'{d}_bounds']) for d in ['R', 'Z', 'T']} if self.is_cyl else {d : asarray(h5f[t0][f'{d}_bounds']) for d in ['X', 'Y', 'Z']}

    if 'T' in self.dims and 'R' in self.dims:
      self.data = vstack((self.data[::-1, :, ::-1, ...], self.data))
      self.abs_error = vstack((self.abs_error[::-1, :, ::-1, ...], self.abs_error))
      self.dims['R'] = concatenate((-1 * self.dims['R'][::-1], self.dims['R']))
      self.dim_bounds['R'] = concatenate((-1 * (self.dim_bounds['R'][1:])[::-1], self.dim_bounds['R']))
  def scale_data(self, scalar):
    self.data *= scalar
    self.abs_error *= scalar
  def set_zero_to_nan(self):
    self.data[self.data == 0] = nan
    self.abs_error[self.data == 0] = nan
  def shift_dims(self, shift_dict):
    for k, v in shift_dict.items():
      if k != 'R' and k != 'T':
        self.dims[k] += v
        self.dim_bounds[k] += v
      else:
        raise ValueError('R and T dimensions cannot be shifted.')
  def scale_dims(self, scalar):
    for k in self.dims:
      if k != 'T':
        self.dims[k] *= scalar
        self.dim_bounds[k] *= scalar
  def slice_by_index(self, slice_dict):
    for dim, slice_val in slice_dict.items():
      if not 0 <= slice_val <= len(self.dims[dim]) - 1:
        raise ValueError(f'Array slice value {slice_val} for dimension {dim} is out of bounds.')
    slices = tuple(slice(v, v+1) if k in self.dims else slice(None) for k, v in slice_dict.items())
    self.data, self.abs_error = squeeze(self.data[slices]), squeeze(self.abs_error[slices])
    self.dims = {k : v for k, v in self.dims.items() if k not in slice_dict}
    self.dim_bounds = {k : v for k, v in self.dim_bounds.items() if k not in slice_dict}
  def slice_by_location(self, slice_dict):
    slices = []
    for k in self.dims:
      if k in slice_dict:
        if self.dim_bounds[k][0] <= slice_dict[k] <= self.dim_bounds[k][-1]:
          idx = int(abs(self.dims[k] - slice_dict[k]).argmin())
          slices.append(slice(idx, idx+1))
        else:
          raise ValueError(f'Specified location {slice_dict[k]} along dimension {k} is out of bounds ({self.dim_bounds[k][0]},{self.dim_bounds[k][-1]})!')
      else:
        slices.append(slice(None))
    slices = tuple(slices)
    self.data, self.abs_error = squeeze(self.data[slices]), squeeze(self.abs_error[slices])
    self.dims = {k : v for k, v in self.dims.items() if k not in slice_dict}
    self.dim_bounds = {k : v for k, v in self.dim_bounds.items() if k not in slice_dict}
  def get_value(self, idx_tuple):
    return ufloat(self.data[idx_tuple], self.abs_error[idx_tuple])
  def sum_along_dim(self, dim):
    axs = list(self.dims.keys()).index(dim)
    self.data = self.data.sum(axis=axs)
    self.abs_error = sqrt(pow(self.abs_error, 2).sum(axis=axs))
    del self.dims[dim]
    del self.dim_bounds[dim]
  def get_bounds(self):
    bounds = []
    for d in self.dim_bounds:
      bounds += [self.dim_bounds[d][0], self.dim_bounds[d][-1]]
    return tuple(bounds)
  def get_abcissa(self):
    return self.dims[list(self.dims.keys())[0]]
  def get_ordinate(self):
    return self.dims[list(self.dims.keys())[1]]

def get_arguments() -> Namespace:
  parser = ArgumentParser(description="A tool to plot mesh tally data")
  plot_group = parser.add_mutually_exclusive_group(required=True)
  plot_group.add_argument(
    '-i', dest='index',
    type=str,
    help='Comma-delimited list of indcies (0-based) at which to slice data. If one index is provided, a 2D plane will be displayed (e.g. x=15). If two indicies are provided, a line will be displayed (e.g. y=5,z=10).',
    metavar=f'D1=<index 1>[,D2=<index 2>]'
  )
  plot_group.add_argument(
    '-v', dest='value',
    type=str,
    help=f'Comma-delimited list of location (in specified units) at which to slice data. If one value is provided, a 2D plane will be displayed. If two values are provided, a line will be displayed.',
    metavar=f'<D1=location 1>[,<D2=location 2>]'
  )
  plot_group.add_argument(
    '-sum', dest='sum_dim',
    help='sum along specified index',
    type=str,
    metavar='<dimension>'
  )
  parser.add_argument(
    '-e', dest='e_index',
    type=int,
    help='Energy index of data to plot. If not provided, data is summed across energy index.',
    metavar='<e index>'
  )
  parser.add_argument(
    '-c', dest='colormap',
    type=str,
    default='nipy_spectral',
    help='Colormap to use for plot. Default: %(default)s',
    metavar='<colormap>'
  )
  value_parser = parser.add_mutually_exclusive_group(required=False)
  value_parser.add_argument(
    '-l', dest='colorbar_label',
    type=str,
    help='Add a colorbar with this label',
    metavar='<colorbar label>'
  )
  value_parser.add_argument(
    '-thr', dest='threshold',
    type=str,
    help='Comma-delimited value:color entries for plotting threshold values. The provided value is the upper bound for that range (with 0 being the implicit lower bound of the lowest range) and the provided color corresponds to that range.',
    metavar='<value1:color1[,value2:color2...]>'
  )
  parser.add_argument(
    '-s', dest='scalar',
    type=float,
    help='Scale data by this number',
    metavar='<scalar>'
  )
  parser.add_argument(
    '-t', dest='title',
    type=str,
    help='Plot title',
    metavar='<title>'
  )
  parser.add_argument(
    '-lb', dest='lower_bound',
    type=float,
    help='Specified lower bound for plot',
    metavar='<lower bound>'
  )
  parser.add_argument(
    '-ub', dest='upper_bound',
    type=float,
    help='Specified upper bound for plot',
    metavar='<upper bound>'
  )
  parser.add_argument(
    '-at', dest='abcissa_tick_interval',
    type=float,
    help='Abcissa major tick interval in specified units',
    metavar='<interval>'
  )
  parser.add_argument(
    '-ot', dest='ordinate_tick_interval',
    type=float,
    help='Ordinate major tick interval in specified units',
    metavar='<interval>'
  )
  parser.add_argument(
    '-al', dest='abcissa_label',
    type=str,
    help='Abcissa label',
    metavar='<label>'
  )
  parser.add_argument(
    '-ol', dest='ordinate_label',
    type=str,
    help='Ordinate label',
    metavar='<label>'
  )
  parser.add_argument(
    '-b', dest='unsampled_color',
    type=str,
    help='Color for unsampled regions',
    metavar='<color>'
  )
  parser.add_argument(
    '-p', dest='plot_data',
    type=str,
    choices=['result', 'aerror', 'rerror'],
    default='result',
    help='Information to plot. Choices=result | aerror (abs error) | rerror (rel error). Default: %(default)s'
  )
  parser.add_argument(
    '-eb', dest='plot_errorbar',
    action='store_true',
    help='Flag to include 1-sigma errorbars on 1D plots'
  )
  parser.add_argument(
    '-pr', dest='probe',
    type=str,
    nargs=1,
    help='Colon-delimited list of coordinates to probe',
    metavar='<pa1,po1[:pa2,po2...]>'
  )
  parser.add_argument(
    '-u', dest='unit',
    type=str,
    choices=['mm', 'cm', 'm', 'km', 'ft'],
    default='cm',
    help='Spatial unit for the plot axes'
  )
  parser.add_argument(
    '-sh', dest='shift',
    type=str,
    help=f'Comma-delimited axis shifts (in specified units) for plot',
    metavar=f'<D1=shift 1>[,<D2=shift 2>]'
  )
  parser.add_argument(
    '--log-scale', dest='logscale',
    action='store_true',
    help='Plot with log-scale'
  )
  parser.add_argument(
    '--no-axis', dest='no_axis',
    action='store_true',
    help='Do not display axis on plot'
  )
  parser.add_argument(
    '--grid', dest='grid',
    type=str,
    nargs='*',
    help='Display grid on plot. Pass comma-delimited keyword arguments to matplotlib grid function.'
  )
  parser.add_argument(
    '-lw', dest='linewidth',
    type=float,
    default=1.0,
    help='Linewidth for 1D plots. Default: %(default)s',
    metavar='<linewidth>'
  )
  parser.add_argument(
    '-hl', dest='hl',
    type=str,
    help='Colon-delimited list of groups of comma-delimited horizontal line information.',
    metavar='<line 1>,[line 2...]'
  )
  parser.add_argument(
    '-vl', dest='vl',
    type=str,
    help='Colon-delimited list of groups of comma-delimited vertical line information.',
    metavar='<line 1>,[line 2...]'
  )
  parser.add_argument(
    '-tx', dest='text',
    type=str,
    help='Colon-delimited list of groups of comma-delimited text annotation information.',
    metavar='<text 1>,[text 2...]'
  )
  parser.add_argument(
    '-ar', dest='arrow',
    type=str,
    help='Colon-delimited list of groups of comma-delimited arrow placement information, starting with four floats designating the start and stop points of the arrow.',
    metavar='<arrow 1>,[arrow 2...]',
    nargs=1
  )
  parser.add_argument(
    'tally',
    type=str,
    help='Number of tally to plot',
    metavar='<tally #>'
  )
  parser.add_argument(
    'output',
    type=str,
    help='Name of output plot file',
    metavar='<plot name>'
  )
  parser.add_argument(
    'files',
    type=str,
    help='Comma-delimited list of meshtal files to process (only supports hdf5). Can plot the sum of tallies. E.g.: 1+2',
    metavar='<file 1>[,<file 2>]'
  )
  return parser.parse_args()

if __name__ == '__main__':
  set_plot_params()
  args, fig, ax = get_arguments(), None, None
  if args.threshold and args.logscale:
    raise ValueError('Log-scale and threshold arguments are incompatible!')
  for fyle in args.files.split(','):
    if ':' in fyle:
      fyle, this_label = fyle.split(':')
    else:
      this_label = None
    with h5file(fyle, 'r') as f:
      mesh_data = MeshDataset(f, args.tally, args.e_index)
      args.tally = args.tally.split('+')[0]
      mesh_data.scale_dims(units[args.unit])
      if args.shift is not None:
        mesh_data.shift_dims({entry.split('=')[0].upper() : float(entry.split('=')[1]) for entry in args.shift.split(',')})
      # Slice the data by index (if index is specified)
      sdict = None
      if args.index is not None:
        sdict = {entry.split('=')[0].upper() : int(entry.split('=')[1]) for entry in args.index.split(',')}
        mesh_data.slice_by_index(sdict)
      # Slice data by position (if location is specfied)
      sdict = None
      if args.value is not None:
        sdict = {entry.split('=')[0].upper() : float(entry.split('=')[1]) for entry in args.value.split(',')}
        mesh_data.slice_by_location(sdict)
      # Sum the data along an axis (if specified)
      if args.sum_dim is not None:
        mesh_data.sum_along_dim(args.sum_dim)
      # Scale the data if requested
      if args.scalar is not None:
        mesh_data.scale_data(args.scalar)
      mesh_data.set_zero_to_nan()
      if len(mesh_data.data.shape) == 2:
        if len(args.files.split(',')) > 1:
          raise ValueError('Error, cannot process multiple meshtal files for 2D plots currently.')
        if not (fig and ax):
          if ('Z' not in mesh_data.dims) and mesh_data.is_cyl:
            fig, ax = subplots(subplot_kw=dict(projection="polar"))
          else:
            fig, ax = subplots()
        # Create the colormap
        if args.threshold:
          cm = ListedColormap([entry.split(':')[1] for entry in args.threshold.split(',')])
          tnorm = BoundaryNorm([0] + [float(entry.split(':')[0]) for entry in args.threshold.split(',')], ncolors=cm.N, clip=True)
        else:
          cm = get_cmap(args.colormap)
          if args.unsampled_color:
            cm.set_bad(color=args.unsampled_color)
        if args.plot_data == 'result':
          pdata = mesh_data.data
        elif args.plot_data == 'aerror':
          pdata = mesh_data.abs_error
        elif args.plot_data == 'rerror':
          pdata = mesh_data.abs_error / mesh_data.data
        # Determine min and max values of the colorbar
        min_data = nanmin(pdata) if args.lower_bound is None else args.lower_bound
        max_data = nanmax(pdata) if args.upper_bound is None else args.upper_bound
        # Change to log-scale if requested
        if args.logscale:
          if ('R' in mesh_data.dims) and ('T' in mesh_data.dims):
            im = ax.pcolormesh(
              2 * pi * mesh_data.dim_bounds['T'],
              mesh_data.dim_bounds['R'][mesh_data.dim_bounds['R'] >= 0.],
              pdata[mesh_data.dims['R'] > 0., ...],
              shading='flat',
              cmap='nipy_spectral',
              norm=LogNorm(vmin=min_data, vmax=max_data)
            )
          else:
            im = ax.imshow(
              pdata.T,
              extent=mesh_data.get_bounds(),
              cmap=cm,
              norm=LogNorm(vmin=min_data, vmax=max_data),
              interpolation='none',
              origin='lower'
            )
        else:
          if ('R' in mesh_data.dims) and ('T' in mesh_data.dims):
            im = ax.pcolormesh(
              2 * pi * mesh_data.dim_bounds['T'],
              mesh_data.dim_bounds['R'][mesh_data.dim_bounds['R'] >= 0.],
              pdata[mesh_data.dims['R'] > 0., ...],
              shading='flat',
              cmap='nipy_spectral',
              vmin=min_data,
              vmax=max_data
            )
          else:
            im = ax.imshow(
              pdata.T,
              extent=mesh_data.get_bounds(),
              cmap=cm,
              vmin=min_data if not args.threshold else None,
              vmax=max_data if not args.threshold else None,
              interpolation='none',
              origin='lower',
              norm=tnorm if args.threshold else None
            )
        # Set aspect ratio to equal so that mesh cells are displayed correctly
        if not (('Z' in mesh_data.dims) and ('T' in mesh_data.dims)):
          ax.set_aspect('equal')
        else:
          ax.set_aspect('auto')
        # Set axis labels appropriately
        if args.no_axis:
          ax.axis('off')
        elif not (('R' in mesh_data.dims and 'T' in mesh_data.dims) and mesh_data.is_cyl):
          if args.abcissa_label is not None:
            ax.set_xlabel(fr'{args.abcissa_label}')
          else:
            ax.set_xlabel(f'Position [{args.unit}]')
          if args.ordinate_label is not None:
            ax.set_ylabel(fr'{args.ordinate_label}')
          else:
            if ('Z' in mesh_data.dims) and ('T' in mesh_data.dims):
              ax.set_ylabel(r'$\theta$ [rev]')
            else:
              ax.set_ylabel(f'Position [{args.unit}]')
        # Add colorbar with specified label if requested
        if args.colorbar_label is not None:
          fig.colorbar(im, label=fr'{args.colorbar_label}', cax=fig.add_axes([ax.get_position().x1+(0.06 if (('Z' not in mesh_data.dims) and mesh_data.is_cyl) else 0.01), ax.get_position().y0, 0.02, ax.get_position().height]))
        # Use scientific notation on colorbar (this is automatic if log-scale)
        if not args.logscale and not (('Z' not in mesh_data.dims) and mesh_data.is_cyl):
          ticklabel_format(useMathText=True)
        if args.probe:
          for loc in args.probe[0].split(':'):
            abcissa, ordinate = mesh_data.get_abcissa(), mesh_data.get_ordinate()
            a, o = loc.split(',')
            a, o = float(a), float(o)
            aidx, oidx = (abs(abcissa - a)).argmin(), (abs(ordinate - o)).argmin()
            ax.scatter([a], [o], marker='D', s=140, color='white', edgecolors='black')
            ax.text(a, o, f'{mesh_data.get_value((aidx, oidx)):.1uS}', ha='center', va='center', color='black', fontsize=4)
      elif len(mesh_data.data.shape) == 1:
        if not (fig and ax):
          fig, ax = subplots()
        abcissa = mesh_data.get_abcissa()
        if args.plot_errorbar:
          ax.errorbar(abcissa, mesh_data.data, yerr=mesh_data.abs_error/2, capsize=args.linewidth, linewidth=args.linewidth, capthick=args.linewidth, label=this_label if this_label else None)
        else:
          ax.plot(abcissa, mesh_data.data, linewidth=args.linewidth, label=this_label if this_label else None)
        ax.set_xlim(mesh_data.get_bounds())
        # Use scientific notation on colorbar (this is automatic if log-scale)
        if args.logscale:
          ax.set_yscale('log')
          ax.yaxis.set_major_locator(LogLocator(base = 10.0))
          ax.yaxis.set_minor_locator(LogLocator(base = 10.0, subs = arange(1.0, 10.0) * 0.1, numticks = 300))
          ax.yaxis.set_minor_formatter(NullFormatter())
        if args.abcissa_label is not None:
          ax.set_xlabel(fr'{args.abcissa_label}')
        else:
          if 'T' in mesh_data.dims:
            ax.set_xlabel(f'Position [rev]')
          else:
            ax.set_xlabel(f'Position [{args.unit}]')
        if args.ordinate_label is not None:
          ax.set_ylabel(fr'{args.ordinate_label}')
        ax.set_ylim((args.lower_bound, args.upper_bound))

  if args.hl:
    args.hl = args.hl.split(':')
    for hl_info in args.hl:
      ax.axhline(**{s.split('=')[0] : float(s.split('=')[1]) if is_float(s.split('=')[1]) else s.split('=')[1] for s in hl_info.split(',')})
  if args.vl:
    args.vl = args.vl.split(':')
    for vl_info in args.vl:
      ax.axvline(**{s.split('=')[0] : float(s.split('=')[1]) if is_float(s.split('=')[1]) else s.split('=')[1] for s in vl_info.split(',')})
  if args.text:
    args.text = args.text.split(':')
    for text_info in args.text:
      ax.text(**{s.split('=')[0] : float(s.split('=')[1]) if is_float(s.split('=')[1]) else s.split('=')[1] for s in text_info.split(',')})
  if args.arrow:
    arrow_opts = args.arrow[0].split(':')
    for arrow_opt in arrow_opts:
      start_x, start_y, stop_x, stop_y, *kwargs = arrow_opt.split(',')
      ax.annotate("", (float(start_x), float(start_y)), (float(stop_x), float(stop_y)), arrowprops={s.split('=')[0] : float(s.split('=')[1]) if is_float(s.split('=')[1]) else s.split('=')[1] for s in kwargs})
  # Set title if requested
  if args.title:
    ax.set_title(fr'{args.title}')
  if args.abcissa_tick_interval:
    ax.xaxis.set_major_locator(MultipleLocator(base=args.abcissa_tick_interval))
  if args.ordinate_tick_interval:
    ax.yaxis.set_major_locator(MultipleLocator(base=args.ordinate_tick_interval))
  # Turn on grid if requested
  if args.grid is not None:
    if len(args.grid) > 0:
      ax.grid(**{s.split('=')[0] : s.split('=')[1] for s in args.grid[0].split(',')})
    else:
      ax.grid(which='both')
  if this_label:
    ax.legend()
  fig.savefig(args.output, bbox_inches='tight', transparent=True)
