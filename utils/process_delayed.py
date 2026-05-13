#!/usr/bin/env python3

from argparse import ArgumentParser
from uncertainties import ufloat
from uncertainties.unumpy import uarray
from numpy import empty, array
import matplotlib.pyplot as plt


class TimeTally:
  def __init__(self, time, data, total, integrate=True):
    self.time, self.total = time, total
    if integrate: # passed time data is assumed to be differential and will be integrated
      tmp_data, tmp_unc = empty(len(data)), empty(len(data)) # temp storage
      self.differential_data = data
      for i in range(len(self.differential_data)):
        tmp_d = (self.differential_data[(len(self.differential_data)-i-1):]).sum()
        tmp_data[len(self.differential_data)-i-1], tmp_unc[len(self.differential_data)-i-1] = tmp_d.nominal_value, tmp_d.std_dev
      self.data = uarray(tmp_data, tmp_unc)
    else:
      self.data = data
  def __sub__(self, other_tally):
    new_tally = TimeTally(self.time, self.data - other_tally.data, self.total - other_tally.total, integrate=False)
    new_tally.differential_data = self.differential_data - other_tally.differential_data
    new_tally.calc_acc()
    return new_tally
  def calc_acc(self):
    tmp_a_data, tmp_a_unc = empty(len(self.data)), empty(len(self.data))
    for i in range(len(self.differential_data)):
      tmp_a = (self.differential_data[:(i+1)]).sum() / self.total
      tmp_a_data[i], tmp_a_unc[i] = tmp_a.nominal_value, tmp_a.std_dev
    self.accumulated_data = uarray(tmp_a_data, tmp_a_unc)

def get_arguments():
  parser = ArgumentParser(description="A tool to compute delayed particle information from MCNP results.")
  parser.add_argument(
    'kcode_results',
    help='MCNP KCODE output file to read',
    type=str,
    metavar='<KCODE output file>'
  )
  parser.add_argument(
    'sdef_delayed',
    help='MCNP SDEF output file (with delayed particles enabled) to read',
    type=str,
    metavar='<SDEF output file w/ delayed>'
  )
  parser.add_argument(
    'sdef_no_delayed',
    help='MCNP SDEF output file (without delayed particles enabled) to read',
    type=str,
    metavar='<SDEF output file w/o delayed>'
  )
  parser.add_argument(
    'tally_number',
    help='MCNP tally number to process (must be the same across all outputs)',
    type=str,
    metavar='<MCNP tally number>'
  )
  parser.add_argument(
    '-p', '--plot',
    help='Name of the plot file.',
    type=str,
    nargs=2,
    metavar='<plot file name> <y-axis label>'
  )
  return parser, parser.parse_args()

def get_time_tally(file_name, tally_number):
  with open(file_name, 'r') as f:
    lines = f.readlines()
    for i, line in enumerate(lines):
      sline = line.split()
      if sline:
        if len(sline) >= 2:
          if sline[0] == '1tally' and sline[1] == tally_number:
            times, vals, vunc = [], [], []
            for j, line2 in enumerate(lines[i+2:]):
              sline2 = line2.split()
              if sline2:
                if sline2[0] == 'time:':
                  [times.append(t) for t in sline2[1:]]
                  next_line = lines[i+2+j+1].split()
                  for k, num in enumerate(next_line):
                    if k % 2 == 0:
                      vals.append(float(num))
                    else:
                      vunc.append(float(num))
              if 'total' in line2:
                break
  times = array([float(t)/6E+09 for t in times[:-1]])
  total, tunc = vals[-1], vunc[-1]
  vals, vunc = array(vals[:-1]), array(vunc[:-1])
  return TimeTally(times, uarray(vals, vals * vunc), ufloat(total, total * tunc))

def get_kcode_tally(file_name, tally_number):
  with open(file_name, 'r') as f:
    lines = f.readlines()
    for i, line in enumerate(lines):
      sline = line.split()
      if sline:
        if len(sline) >= 2:
          if sline[0] == '1tally' and sline[1] == tally_number:
            tline = lines[i+12].split()
            if tline[0] == 'cell':
              tline = lines[i+13].split()
            v, u = float(tline[0]), float(tline[1])
            return ufloat(v, v * u)

if __name__ == '__main__':
  argp, args = get_arguments()
  sdef_delayed, sdef_prompt = (get_time_tally(f, args.tally_number) for f in [args.sdef_delayed, args.sdef_no_delayed])
  kcode = get_kcode_tally(args.kcode_results, args.tally_number)
  sdef_delayed_only = sdef_delayed - sdef_prompt
  delayed_fraction = (sdef_delayed.total - sdef_prompt.total) / sdef_prompt.total
  print(f'Delayed to prompt ratio: {delayed_fraction:.1uS}')
  if args.plot:
    fig, ax = plt.subplots()
    plot_dat = sdef_delayed_only.accumulated_data * delayed_fraction * kcode
    plot_vals = array([v.nominal_value for v in plot_dat])
    ax.plot(sdef_delayed_only.time, plot_vals)
    ax.set_xlabel('Time (minutes)')
    ax.set_ylabel(fr'{args.plot[1]}')
    ax.set_xlim([0, max(sdef_delayed_only.time)])
    ax.set_ylim([0, max(plot_vals)])
    ax.grid()
    plt.savefig(args.plot[0], bbox_inches='tight')
