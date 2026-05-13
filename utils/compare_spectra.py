#! /usr/bin/env python3

from argparse import ArgumentParser
from numpy import asarray
from math import isclose
from uncertainties import ufloat
from uncertainties.unumpy import uarray, pow
from matplotlib import ticker
from matplotlib.pyplot import subplots, hist, savefig, errorbar
from typing import Union
# Local modules
from mcnp_utilities.utils.extract_spectrum import SpectrumTally
from mcnp_utilities.utils.plot_spectra import read_tally_file


def compare_spectra(spec_list: list[SpectrumTally], plot_files: Union[list[str], None] = None, plot_labels: Union[list[str], None] = None, absolute: bool = False, xlims: Union[list[float], None] = None, ylims: Union[list[float], None] = None) -> None:
  ref_spec = spec_list[0]
  ref_spec.tally_data /= sum(ref_spec.tally_data) # Normalize reference
  for idx, spec in enumerate(spec_list):
    if idx != 0:
      spec.tally_data /= sum(spec.tally_data) # Normalize comparison data
      print(f'Comparing spectrum at index {idx} to reference at index 0...')
      assert spec.bin_number == ref_spec.bin_number
      assert all(isclose(x, y) for x, y in zip(spec.energy_bin_edges, ref_spec.energy_bin_edges))
      rel_diffs, rel_diffs_unc = [], []
      for e_idx, (comp, ref_d) in enumerate(zip(spec.tally_data, ref_spec.tally_data)):
        rel_diff = (100. * abs(comp - ref_d) / ref_d if not isclose(ref_d.nominal_value, 0.) else ufloat(0., 0.)) if absolute else (100. * (comp - ref_d) / ref_d if not isclose(ref_d.nominal_value, 0.) else ufloat(0., 0.))
        rel_diffs.append(rel_diff.nominal_value)
        rel_diffs_unc.append(rel_diff.std_dev)
        print(f'  - rel diff for index {e_idx} : {rel_diff:.1uS}%')
      diffs_array = uarray(rel_diffs, rel_diffs_unc)
      rms_val = pow((diffs_array**2).mean(), 0.5)
      print(f'avg relative difference: {(diffs_array).mean():.1uS}%')
      print(f'RMS of relative difference: {rms_val:.1uS}%')
      if plot_files:
        dat, err = asarray([d.nominal_value for d in diffs_array]), asarray([d.std_dev for d in diffs_array])
        fig, ax = subplots()
        (_,_,patches) = hist(
          x=ref_spec.bin_midpoints,
          bins=ref_spec.energy_bin_edges,
          weights=dat,
          linewidth=0.5,
          histtype='step',
          linestyle='solid',
          color='r'
        )

        (_,caps,_) = errorbar(
          ref_spec.bin_midpoints,
          dat,
          yerr=err/2,
          fmt='',
          ls='none',
          linewidth=0.5,
          color='r'
        )

        ax.grid()
        ax.set_xscale('log')
        ax.set_xlim(tuple(xlims)) if xlims is not None else ax.set_xlim((ref_spec.energy_bin_edges[0], ref_spec.energy_bin_edges[-1]))
        if ylims is not None:
          ax.set_ylim(tuple(ylims))
        ax.set_xlabel('Energy [MeV]')
        if plot_labels:
          ax.set_ylabel(fr'{plot_labels[idx-1]}')
        else:
          ax.set_ylabel('Absolute Relative Difference [%]' if absolute else 'Relative Difference [%]')

        # Add all decade ticks and sub-decade ticks
        locmaj = ticker.LogLocator(base=10,numticks=300)
        ax.xaxis.set_major_locator(locmaj)
        locmin = ticker.LogLocator(base=10.0,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=300)
        ax.xaxis.set_minor_locator(locmin)
        ax.xaxis.set_minor_formatter(ticker.NullFormatter())

        savefig(plot_files[idx-1], bbox_inches='tight')


def get_arguments():
  parser = ArgumentParser(description='A tool to compare energy spectra from MCNP output files')
  parser.add_argument(
    '-p', dest='plots',
    type=str,
    nargs='+',
    help=r'Name(s) of plot file(s) for plotting relative error (percent) as a function of energy',
    metavar='<plot name>'
  )
  parser.add_argument(
    '-l', dest='labels',
    type=str,
    nargs='+',
    help='y-axis labels for plots',
    metavar='<y-axis label>'
  )
  parser.add_argument(
    '-a', dest='absolute',
    action='store_true',
    default=False,
    help='Flag to plot absolute differences. Default: %(default)s'
  )
  parser.add_argument(
    '-x', dest='xlims',
    type=float,
    nargs=2,
    help='abcissa lower and upper bounds',
    metavar=('<lower bound>', '<upper bound>')
  )
  parser.add_argument(
    '-y', dest='ylims',
    type=float,
    nargs=2,
    help='ordinate lower and upper bounds',
    metavar=('<lower bound>', '<upper bound>')
  )
  parser.add_argument(
    'files',
    type=str,
    nargs='+',
    help='Spectrum files to read (first entry is treated as reference)',
    metavar='<spectrum file>'
  )
  return parser, parser.parse_args()

if __name__ == '__main__':
  argp, args = get_arguments()
  if args.plots:
    if len(args.plots) != (len(args.files) - 1):
      argp.error(f'Number of plot file names {(len(args.plots))} is not equal to the number of comparisons ({len(args.files) - 1})!')
  if args.labels and not args.plots:
    argp.error('Cannot specify plot labels without plot option!')
  if args.plots and args.labels:
    if len(args.plots) != len(args.labels):
      argp.error(f'Number of plot file names {(len(args.plots))} is not equal to the number of plot labels ({len(args.labels)})!')

  compare_spectra([read_tally_file(s) for s in args.files], args.plots, args.labels, args.absolute, args.xlims, args.ylims)
