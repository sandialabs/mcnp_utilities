#! /usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib import colors
from argparse import ArgumentParser, Action
from pandas import DataFrame, ExcelWriter
from numpy import argmin, argmax
from uncertainties import ufloat
from h5py import File as h5file
from numpy import asarray
# Local modules
from mcnp_utilities.utils.get_rx import get_rho
from mcnp_utilities.lib.basic_tools import hex_2_cart, set_plot_params


set_plot_params()

class CommaSplitAction(Action):
  def __call__(self, parser, namespace, values, option_string=None):
    setattr(namespace, self.dest, values.split(","))

class ColonandCommaSplitAction(Action):
  def __call__(self, parser, namespace, values, option_string=None):
    setattr(namespace, self.dest, [v.split(",") for v in values.split(":")])

class PlusandColonandCommaSplitAction(Action):
  def __call__(self, parser, namespace, values, option_string=None):
    setattr(namespace, self.dest, [[v.split(",") for v in v2.split(":")] for v2 in values.split('+')])

class RodTally:
  def __init__(self, val, px, py):
    self.value, self.xpos, self.ypos, self.pf = val, px, py, 0.

class Core:
  def __init__(self, rod_dictionary, name):
    self.rod_dict, self.power, self.num_rods, self.avg_power, self.name = rod_dictionary, 0., 0, 0., name
  def compute_power(self):
    tmp_dict, self.power, self.num_rods = {}, 0., 0
    for k in self.rod_dict:
      tmp_dict[k] = [self.rod_dict[k][x].value for x in range(len(self.rod_dict[k]))]
      self.power += sum(tmp_dict[k])
      self.num_rods += len(self.rod_dict[k])
    del tmp_dict
  def compute_pf(self, total_power, total_rods):
    avg_power = total_power / total_rods
    for k in self.rod_dict:
      for rod in self.rod_dict[k]:
        rod.pf = rod.value / avg_power
  def print_pf(self):
    max_of_max = ufloat(0., 0.)
    print(f'\nPF data for {self.name}:')
    for k in self.rod_dict:
      print(f"\n{k}\n{'-'*len(k)}")
      print(f"{'x':^5}|{'y':^5}|{'PF':^6}")
      max_pf = ufloat(0., 0.)
      for rod in self.rod_dict[k]:
        print(f"{rod.xpos:^5n}|{rod.ypos:^5n}|{rod.pf:^6.1uS}")
        if rod.pf > max_pf:
          max_pf = rod.pf
        if rod.pf > max_of_max:
          max_of_max = rod.pf
      print(f"Max PF for {k}: {max_pf:.1uS}")
    print(f"\nMax PF for {self.name}: {max_of_max:.1uS}")
    self.avg_power = self.power / self.num_rods
    print(f"Average Power for {self.name}: {self.avg_power:.1uS}")

def get_plt_ratio(core_list, xtr, ytr):
  x_vals, y_vals, max_x, max_y = [], [], 0, 0
  for core, x, y in zip(core_list, xtr, ytr):
    x_tr, y_tr = float(x), float(y)
    for k in core.rod_dict:
      for rod in core.rod_dict[k]:
        max_x = max(max_x, rod.xpos)
        max_y = max(max_y, rod.ypos)
        x_loc, y_loc = hex_2_cart(rod.xpos, rod.ypos, core.pitch)
        x_loc += x_tr
        y_loc += y_tr
        x_vals.append(x_loc)
        y_vals.append(y_loc)
    min_xl, max_xl = min(x_vals), max(x_vals)
    min_yl, max_yl = min(y_vals), max(y_vals)
    tmpx = abs(min_xl) + abs(max_xl) + 1
    tmpy = abs(min_yl) + abs(max_yl) + 1
  return 68. / max(tmpx, tmpy)

def plot_2d_hex(core_dict, xtr, ytr, percent, lbound, ubound, font_scalar, marker_scalar, nominal_only=None, highlight=None, no_number=False):
  fig, axes = plt.subplots(nrows=1, ncols=len(core_dict))
  axes = [axes] if len(core_dict) == 1 else axes
  plot_data = {}
  plot_center = 0 if percent else 1
  x_vals, y_vals = {}, {}
  for i, (ax, (f, core_list)) in enumerate(zip(axes, core_dict.items())):
    fs = float(font_scalar[i])
    ax.set_aspect('equal')
    plt_ratio = get_plt_ratio(core_list, xtr[i], ytr[i])
    x_vals[f], y_vals[f], pf_vals = [], [], []
    for core, x, y in zip(core_list, xtr[i], ytr[i]):
      x_tr, y_tr = float(x), float(y)
      for k in core.rod_dict:
        for rod in core.rod_dict[k]:
          x_loc, y_loc = hex_2_cart(rod.xpos, rod.ypos, core.pitch)
          x_loc += x_tr
          y_loc += y_tr
          x_vals[f].append(x_loc)
          y_vals[f].append(y_loc)
          pf_vals.append(rod.pf)
    minl, maxl, min_pf, max_pf, pf_percent = argmin(pf_vals), argmax(pf_vals), min(pf_vals), max(pf_vals), []
    if not no_number:
      for i, (x, y, pf) in enumerate(zip(x_vals[f], y_vals[f], pf_vals)):
        if highlight is not None:
          if pf > highlight:
            clr = 'yellow'
          else:
            if (0.75 * (1. - min_pf) + min_pf) < pf < (0.25 * (max_pf - 1.) + 1.):
              clr = 'black'
            else:
              clr = 'white'
        else:
          if i == maxl:
            clr = 'magenta'
          elif i == minl:
            clr = 'cyan'
          elif (0.75 * (1. - min_pf) + min_pf) < pf < (0.25 * (max_pf - 1.) + 1.):
            clr = 'black'
          else:
            clr = 'white'
        if percent:
          pf_percent.append(100 * (pf - 1))
          if nominal_only:
            ax.text(x, y - 0.08 * core.pitch, f"{pf_percent[-1].nominal_value:.{nominal_only}f}%", fontsize=fs * round(4 * plt_ratio), ha='center', c=clr, weight='bold')
          else:
            ax.text(x, y - 0.05 * core.pitch, f"{pf_percent[-1]:.1uS}%", fontsize=fs * round(2 * plt_ratio), ha='center', c=clr, weight='bold')
        else:
          if nominal_only:
            ax.text(x, y - 0.11 * core.pitch, f"{pf.nominal_value:.{nominal_only}f}", fontsize=fs * round(4 * plt_ratio), ha='center', c=clr, weight='bold')
          else:
            ax.text(x, y - 0.05 * core.pitch, f"{pf:.1uS}", fontsize=fs * round(2 * plt_ratio), ha='center', c=clr, weight='bold')
    plot_data[f] = [p.nominal_value for p in pf_percent] if percent else [p.nominal_value for p in pf_vals]
  l_bound = lbound if lbound is not None else min([min(d) for d in plot_data.values()])
  u_bound = ubound if ubound is not None else max([max(d) for d in plot_data.values()])
  divnorm=colors.TwoSlopeNorm(vmin=l_bound, vcenter=plot_center, vmax=u_bound)
  for i, (ax, f) in enumerate(zip(axes, args.file)):
    sc = ax.scatter(x_vals[f], y_vals[f], c=plot_data[f], cmap='seismic', norm=divnorm, marker='h', s=float(marker_scalar[i]) * 285. * plt_ratio**2)
    ax.set_frame_on(False)
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
  for ax in axes:
    ax.set_xlim((min([ax2.get_xlim()[0] for ax2 in axes]), max([ax2.get_xlim()[1] for ax2 in axes])))
    ax.set_ylim((min([ax2.get_ylim()[0] for ax2 in axes]), max([ax2.get_ylim()[1] for ax2 in axes])))
  cbar = plt.colorbar(
    sc,
    orientation="vertical",
    label='% Deviation from Avg. Power' if percent else 'Peaking Factor',
    cax=fig.add_axes(
      [
        axes[-1].get_position().x1+0.01,
        axes[-1].get_position().y0,
        0.02,
        axes[-1].get_position().height
      ]
    )
  )
  if not no_number:
    cbar.set_ticks([l_bound, plot_center, u_bound])
    if percent:
      cbar.set_ticklabels([f"{l_bound:.0f}", f"{plot_center:.0f}", f"{u_bound:.0f}"])
    else:
      cbar.set_ticklabels([f"{l_bound:.2f}", f"{plot_center:.2f}", f"{u_bound:.2f}"])
  else:
    cbar.set_ticks([])
  for i, (ax, f) in enumerate(zip(axes, args.file)):
    if args.reactivity:
      rho = get_rho(f)
      ax.set_title(rf"$\rho$ = {rho:.1uS}\$".replace('-', u'\u2013'))
    if args.sub_label:
      ax.axes.get_xaxis().set_visible(True)
      ax.tick_params(axis='x', bottom=False)
      ax.axes.get_xaxis().set_ticks([])
      ax.set_xlabel(args.sub_label[i])
    ax.autoscale_view()
  fig.subplots_adjust(wspace=0, hspace=0)
  if args.map.endswith('.png'):
    plt.savefig(args.map, bbox_inches='tight', dpi=500, transparent=True)
  else:
    plt.savefig(args.map, bbox_inches='tight', transparent=True)

def write_excel(cores, name):
  with ExcelWriter(name) as writer:
    for core in cores:
      for k in core.rod_dict:
        ddict = {k : [] for k in ['x', 'y', 'PF']}
        for rod in core.rod_dict[k]:
          ddict['x'].append(rod.xpos)
          ddict['y'].append(rod.ypos)
          ddict['PF'].append(f'{rod.pf:.1uS}')
        DataFrame(ddict).to_excel(writer, sheet_name=k, index=False)

def write_hdf5(cores, name):
  with h5file(name, 'w') as f:
    for core in cores:
      for k in core.rod_dict:
        g = f.create_group(k)
        x_pos, y_pos, pf_data, pf_err = [], [], [], []
        for row, rod in enumerate(core.rod_dict[k]):
          x_pos.append(rod.xpos)
          y_pos.append(rod.ypos)
          pf_data.append(rod.pf.nominal_value)
          pf_err.append(rod.pf.std_dev)
        g['x_coordinate'] = asarray(x_pos)
        g['y_coordinate'] = asarray(y_pos)
        g['peaking_factor'] = asarray(pf_data)
        # g['peaking_factor_stddev'] = asarray(pf_err)

def read_rod_data(inp_file, fuel_dictionary, tally_number):
  rod_dictionary, in_tally = {}, False
  for v in fuel_dictionary.values():
    rod_dictionary[v] = []
  with open(inp_file, 'r', errors='ignore') as f:
    for i, line in enumerate(f):
      sline = line.split()
      if len(sline) > 2:
        if sline[0] == "1tally" and sline[1] == tally_number and "print table 30" not in line:
          in_tally = True
        if line.startswith(" cell (") and in_tally:
          for k in fuel_dictionary:
            if line.startswith(f" cell ({k}<"):
              lline = line.split("[")[1].split("]")[0].split()
              x, y = int(lline[0]), int(lline[1])
              for j in range(2):
                line = next(f)
              tline = line.split()
              v, e = float(tline[0]), float(tline[1])
              v = ufloat(v, v*e)
              if v > 0.:
                rod_dictionary[fuel_dictionary[k]].append(RodTally(v, x, y))
        if in_tally and "cell union total" in line:
          break
  return rod_dictionary

def getArguments():
  parser = ArgumentParser(description="A tool to compute core peaking factors using an MCNP output file.")
  parser.add_argument(
    '-n', dest='names',
    help='Cell names, comma-delimited, with tallies delimited by a colon, and files delimited by a +',
    type=str,
    action=PlusandColonandCommaSplitAction,
    metavar='<name1[,name2...]>'
  )
  parser.add_argument(
    '-c', dest='core_names',
    help='Comma-delimited list of core names, with a colon delimiting files',
    type=str,
    action=ColonandCommaSplitAction,
    metavar='<corename1[,corename2...]>'
  )
  parser.add_argument(
    '-m', dest='map',
    help='Name of peaking factor plot file.',
    type=str,
    metavar='<plot file name>'
  )
  parser.add_argument(
    '-e', dest='excel',
    help='Prefix of Excel file name to which to write (".xlsx" is appended to given name). Must be comma-delimited list if mulitple files are specified.',
    type=str,
    action=CommaSplitAction,
    metavar='<file name prefix 1[,file name prefix 2 ...]>'
  )
  parser.add_argument(
    '-h5', dest='hdf5',
    help='Prefix of HDF5 file name to which to write (".h5" is appended to given name). Must be comma-delimited list if mulitple files are specified.',
    type=str,
    action=CommaSplitAction,
    metavar='<file name prefix 1[,file name prefix 2 ...]>'
  )
  parser.add_argument(
    '-f', dest='print_power_fractions',
    help='Option to print power fractions of each core for a given file (one core is associated with one tally, and fractions are printed in the same order as tallies)',
    action='store_true'
  )
  parser.add_argument(
    '-%', dest='percent',
    help='Flag to use percent deviation from average core power on plot. (default: %(default)s)',
    action='store_true',
    default=False
  )
  parser.add_argument(
    '-nn', dest='no_number',
    help='Flag to omit numbers on the plot and the colorbar. (default: %(default)s)',
    action='store_true',
    default=False
  )
  parser.add_argument(
    '-x', dest='xtranslation',
    help='Comma-delimited list of x-axis core translations, with a colon delimiting files',
    type=str,
    action=ColonandCommaSplitAction,
    metavar='<tr1[,tr2...]>'
  )
  parser.add_argument(
    '-y', dest='ytranslation',
    help='Comma-delimited list of y-axis core translations, with a colon delimiting files',
    type=str,
    action=ColonandCommaSplitAction,
    metavar='<tr1[,tr2...]>'
  )
  parser.add_argument(
    '-p', dest='pitch',
    help='Comma-delimited list of core pitches, with files delimited by a colon (optional but HIGHLY recommended for PF maps).',
    type=str,
    action=ColonandCommaSplitAction,
    metavar='<tr1[,tr2...]>'
  )
  parser.add_argument(
    '-l', dest='lower_bound',
    help='Lower bound for PF map.',
    type=float,
    metavar='<lower bound>'
  )
  parser.add_argument(
    '-u', dest='upper_bound',
    help='Upper bound for PF map.',
    type=float,
    metavar='<upper bound>'
  )
  parser.add_argument(
    '-hl', dest='highlight',
    help='Highlight elements above this value on the plot',
    type=float,
    metavar='<highlight value>'
  )
  parser.add_argument(
    '-r', dest='reactivity',
    help='Add reactivity of model to PF map (can only be used with map option). Default: %(default)s',
    action='store_true'
  )
  parser.add_argument(
    '-s', dest='sub_label',
    help='Add a label below each core (comma-delimited by file)',
    action=CommaSplitAction,
    metavar='<label 1[,label 2 ...]>',
    type=str
  )
  parser.add_argument(
    '-no', dest='nominal_only',
    help='If provided, uncertainties of PFs are omitted on PF map in each element location. This value is the number of decimal places used when printing PF nominal values.',
    type=int,
    metavar='<# decimal places>'
  )
  parser.add_argument(
    '-fs', dest='font_scalar',
    help='Value by which to scale font size of element PFs on PF maps (comma-delimited by file). Default: %(default)s',
    type=str,
    action=CommaSplitAction,
    metavar='<scalar 1[,scalar 2 ...]>'
  )
  parser.add_argument(
    '-ms', dest='marker_scalar',
    help='Value by which to scale marker of element PFs on PF maps (comma-delimited by file). Default: %(default)s',
    type=str,
    action=CommaSplitAction,
    metavar='<scalar 1[,scalar 2 ...]>'
  )
  parser.add_argument(
    'file',
    help='MCNP output file(s) to process (comma-delimited)',
    type=str,
    action=CommaSplitAction,
    metavar='<file 1>[,file 2 ...]'
  )
  parser.add_argument(
    'tallies',
    help='Comma-delimited list of tallies to process, with files delimited by a colon',
    type=str,
    action=ColonandCommaSplitAction,
    metavar='<tally1[,tally2 ...]>'
  )
  parser.add_argument(
    'cells',
    help='Cell numbers (comma-delimited) with tallies delimited by a colon and files delimited by a +',
    type=str,
    action=PlusandColonandCommaSplitAction,
    metavar='<cell1[,cell2 ...]>'
  )
  return parser, parser.parse_args()

def error_check(argp, argo):
  n_tally_lists, n_cell_lists, n_files = len(argo.tallies), len(args.cells), len(argo.file)
  for o, s in zip([n_cell_lists, n_tally_lists], ['cell', 'tally']):
    if o != n_files:
      argp.error(f'Number of file {s} lists ({o}) is not equal to the number of files ({n_files})!')
  if n_cell_lists == n_files:
    for i, tally_cell_list in enumerate(argo.cells):
      if len(tally_cell_list) != len(argo.tallies[i]):
        argp.error(f'Number of file cell lists ({len(tally_cell_list)}) for file index {i} is not equal to the number of tally lists ({len(argo.tallies[i])})!')
  if args.sub_label:
    if n_files != len(args.sub_label):
      argp.error(f'Number of sublabels ({len(args.sub_label)}) is not equal to the number of files ({n_files})!')
  if argo.names:
    if len(argo.names) != n_files:
      argp.error(f'Number of file name lists ({len(argo.names)}) is not equal to the number of files ({n_files})!')
    for i, tally_name_list in enumerate(argo.names):
      if len(tally_name_list) != len(argo.tallies[i]):
        argp.error(f'Number of file name lists ({len(tally_name_list)}) for file index {i} is not equal to the number of tallies ({len(argo.tallies[i])})!')
      for j, cell_name_list in enumerate(argo.names[i]):
        if len(cell_name_list) != len(argo.cells[i][j]):
          argp.error(f'Number of cell names ({len(cell_name_list)}) for file index {i} and tally index {j} is not equal to the number of cells ({len(argo.cells[i][j])})!')
  if argo.core_names:
    if len(argo.core_names) != n_files:
      argp.error(f'Number of file core name lists ({len(argo.core_names)}) is not equal to the number of files ({n_files})!')
    for i, core_name_list in enumerate(argo.core_names):
      if len(core_name_list) != len(argo.tallies[i]):
        argp.error(f'Number of core names ({len(core_name_list)}) for file index {i} is not equal to the number of tallies ({len(argo.tallies[i])})!')
  for l, s in zip([argo.font_scalar, argo.marker_scalar], ['font scalar', 'marker scalar']):
    if l:
      if len(l) != n_files:
        argp.error(f'Number of {s} entries ({len(l)}) is not equal to the number of files ({n_files})!')
  f_opts = [argo.pitch, argo.xtranslation, argo.ytranslation]
  tmp_strs = ['pitch', 'x-translation', 'y-translation']
  for o, s in zip(f_opts, tmp_strs):
    if o:
      if len(o) != n_files:
        argp.error(f'Number of file {s} lists ({len(o)}) is not equal to the number of files ({n_files})!')
      for i, l in enumerate(o):
        if len(l) != len(argo.tallies[i]):
          argp.error(f'Number of core {s}s ({len(l)}) for file index {i} is not equal to the number of tallies ({len(argo.tallies[i])})!')
  if not argo.map:
    map_opts = [argo.percent, argo.reactivity, argo.xtranslation, argo.ytranslation, (argo.lower_bound or argo.upper_bound), argo.font_scalar, argo.marker_scalar, argo.pitch, argo.nominal_only]
    tmp_strs = ['percent option (-%)', 'reactivity option (-r)', 'x-translation option (-x)', 'y-translation option (-y)', 'bound options (-l and -u)', 'font scalar option (-fs)', 'marker scalar option (-ms)', 'pitch option (-p)', 'nominal only option (-no)']
    for o, s in zip(map_opts, tmp_strs):
      if o:
        argp.error(f'Cannot use {s} without map option (-m)!')
    if args.highlight is not None:
      argp.error(f'Cannot use highlight option (-hl) without map option (-m)!')
  for o, s in zip([argo.excel, argo.hdf5], ['excel', 'HDF5']):
    if o:
      if len(o) != n_files:
        argp.error(f'Number of {s} file prefixes ({len(o)}) is not equal to the number of files ({n_files})!')

if __name__ == "__main__":
  aparser, args = getArguments()
  error_check(aparser, args)
  args.names = {f : [[f'{t}-{c}' for c in args.cells[j][i]] for i, t in enumerate(args.tallies[j])] for j, f in enumerate(args.file)} if args.names is None else {f : [args.names[j][i] for i, _ in enumerate(args.tallies[j])] for j, f in enumerate(args.file)}
  fuel_dicts = {f : [{k : v for k, v in zip(c, n)} for c, n in zip(args.cells[i], args.names[f])] for i, f in enumerate(args.file)}

  args.xtranslation = [[0. for n in args.tallies[i]] for i in range(len(args.file))] if args.xtranslation is None else args.xtranslation
  args.ytranslation = [[0. for n in args.tallies[i]] for i in range(len(args.file))] if args.ytranslation is None else args.ytranslation

  rod_dicts = {f : [read_rod_data(f, d, n) for d, n in zip(fuel_dicts[f], args.tallies[i])] for i, f in enumerate(args.file)}
  cores = {f : [Core(rd, args.core_names[j][i]) for i, rd in enumerate(rod_dicts[f])] for j, f in enumerate(args.file)} if args.core_names else {f : [Core(rd, args.tallies[j][i]) for i, rd in enumerate(rod_dicts[f])] for j, f in enumerate(args.file)}
  if args.pitch:
    for (f, c), pitch_list in zip(cores.items(), args.pitch):
      for core, p in zip(c, pitch_list):
        core.pitch = float(p)
  else:
    for f, c in cores.items():
      for core in c:
        core.pitch = 1.
  [[core.compute_power() for core in c] for c in cores.values()]
  tot_power = {f : sum([core.power for core in cores[f]]) for f in args.file}
  tot_rods = {f : sum([core.num_rods for core in cores[f]]) for f in args.file}
  [[core.compute_pf(tot_power[f], tot_rods[f]) for core in c] for f, c in cores.items()]
  for f in args.file:
    print(f'\nPeaking Factor Data for {f}:')
    [core.print_pf() for core in cores[f]]

  if args.map:
    plot_2d_hex(cores,
      args.xtranslation,
      args.ytranslation,
      True if args.percent else False,
      args.lower_bound if args.lower_bound else None,
      args.upper_bound if args.upper_bound else None,
      args.font_scalar if args.font_scalar else [1.0 for i in range(len(cores))],
      args.marker_scalar if args.marker_scalar else [1.0 for i in range(len(cores))],
      args.nominal_only,
      args.highlight if args.highlight is not None else None,
      args.no_number
    )

  if args.print_power_fractions:
    for j, f in enumerate(args.file):
      print(f'\nPower generation fractions for {f}:')
      for i, core in zip(args.tallies[j], cores[f]):
        p_frac = core.power / tot_power[f]
        print(f'Tally {i} : {p_frac:.1uS}')
  if args.excel:
    [write_excel(cores[f], f'{e}.xlsx') for f, e in zip(args.file, args.excel)]
  if args.hdf5:
    [write_hdf5(cores[f], f'{h}.h5') for f, h in zip(args.file, args.hdf5)]
