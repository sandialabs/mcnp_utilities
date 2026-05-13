#! /usr/bin/env python3

from argparse import ArgumentParser
from numpy import sum
from matplotlib import ticker
from matplotlib.pyplot import subplots, hist, errorbar, rcParams
from pickle import load
# Local modules
from mcnp_utilities.utils.extract_spectrum import SpectrumTally
from mcnp_utilities.lib.basic_tools import set_plot_params


def get_arguments():
  parser = ArgumentParser(description='A tool to plot fluence energy spectra')
  parser.add_argument(
    'output_name',
    type=str,
    help='Name of the output plot file',
    metavar='<plot name>'
  )
  parser.add_argument(
    'tally_files',
    help='Tally files to plot',
    nargs='+',
    metavar='<file>'
  )
  parser.add_argument(
    '-l', dest='labels',
    nargs='+',
    type=str,
    help='Plot labels',
    metavar='<plot label>'
  )
  parser.add_argument(
    '-c', dest='colors',
    nargs='+',
    type=str,
    help='Plot colors',
    metavar='<plot color string>'
  )
  parser.add_argument(
    '-s', dest='styles',
    nargs='+',
    type=str,
    help='Plot linestyles',
    metavar='<plot linestyle string>'
  )
  parser.add_argument(
    '-y', dest='y_label',
    type=str,
    help='y-axis label for plot. Default: %(default)s',
    default=r'Fluence $\left[\frac{\mathrm{\#}}{\mathrm{cm}^2 \cdot \mathrm{MJ}}\right]$ per lethargy',
    metavar='<plot label>'
  )
  parser.add_argument(
    '-p', dest='legend_position',
    type=str,
    help='Legend position string. Default: %(default)s',
    default='best',
    metavar='<plot legend position string>'
  )
  parser.add_argument(
    '-n', dest='normalize',
    action='store_true',
    help='Flag indicating whether to normalize tally data to unity (removes ordinate tick labels from plots)'
  )
  parser.add_argument(
    '-e', dest='errorbars',
    action='store_true',
    help='Flag indicating whether to include errorbars on plot'
  )
  parser.add_argument(
    '-a', dest='abcissa',
    type=float,
    nargs=2,
    help='abcissa lower and upper bounds',
    metavar=('<lower bound>', '<upper bound>')
  )
  parser.add_argument(
    '-o', dest='ordinate',
    type=float,
    nargs=2,
    help='ordinate lower and upper bounds',
    metavar=('<lower bound>', '<upper bound>')
  )
  return parser, parser.parse_args()

def read_tally_file(file_name: str) -> SpectrumTally:
  with open(file_name, 'rb') as f:
    return load(f)

if __name__ == '__main__':
  set_plot_params()
  argp, args = get_arguments()
  if args.colors:
    if args.labels:
      if len(args.labels) != len(args.tally_files):
        argp.error('# of label strings must be the same as the number of tally files!')
    if len(args.colors) != len(args.tally_files):
      argp.error('# of color strings must be the same as the number of tally files!')
    if args.styles:
      if len(args.styles) != len(args.tally_files):
        argp.error('# of linestyle strings must be the same as the number of tally files!')
  fig, ax = subplots()
  cc = rcParams['axes.prop_cycle'].by_key()['color']
  for i, f in enumerate(args.tally_files):
    tal = read_tally_file(f)
    results = tal.tally_data

    if args.normalize:
      rawtally = results * tal.u             # undo per lethargy first
      rawsum   = sum(rawtally)               # total raw tally sum value, not total per lethargy sum value
      results  = rawtally / (rawsum * tal.u) # normalized per lethargy tally

    (_,_,patches) = hist(
      x=tal.bin_midpoints,
      bins=tal.energy_bin_edges,
      weights=[v.nominal_value for v in results],
      linewidth=0.5,
      histtype='step',
      linestyle=args.styles[i] if args.styles else 'solid',
      label=fr'{args.labels[i]}' if args.labels else None,
      color=args.colors[i] if args.colors else cc[i]
    )
    if args.errorbars:
      (_,caps,_) = errorbar(
        tal.bin_midpoints,
        [v.nominal_value for v in results],
        yerr=[e.std_dev for e in results],
        fmt='',
        ls='none',
        linewidth=0.5,
        color=args.colors[i] if args.colors else cc[i]
      )
      for cap in caps:
        cap.set_color(args.colors[i] if args.colors else cc[i])
  ax.set_xscale('log')
  ax.set_xlabel('Energy [MeV]')
  if args.abcissa is not None:
    ax.set_xlim(tuple(args.abcissa))
  if args.ordinate is not None:
    ax.set_ylim(tuple(args.ordinate))
  if args.normalize:
    ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
  ax.set_ylabel(fr'{args.y_label}')
  ax.set_axisbelow(True)
  ax.grid()

  # Add all decade ticks and sub-decade ticks
  locmaj = ticker.LogLocator(base=10,numticks=300)
  ax.xaxis.set_major_locator(locmaj)
  locmin = ticker.LogLocator(base=10.0,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=300)
  ax.xaxis.set_minor_locator(locmin)
  ax.xaxis.set_minor_formatter(ticker.NullFormatter())

  if args.labels:
    ax.legend(loc=args.legend_position)

  ax.ticklabel_format(axis='y', useMathText=True)
  fig.savefig(args.output_name, bbox_inches='tight')
