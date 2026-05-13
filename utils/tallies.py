#! /usr/bin/env python3

from uncertainties import ufloat
from argparse import ArgumentParser, Action
from collections import defaultdict
# Local modules
from mcnp_utilities.utils.get_rx import get_keff


class NestedDefaultDict(defaultdict):
  def __init__(self, *args, **kwargs):
    super(NestedDefaultDict, self).__init__(NestedDefaultDict, *args, **kwargs)
  def __repr__(self):
    return repr(dict(self))

class ColonSplitAction(Action):
  def __call__(self, parser, namespace, values, option_string=None):
    setattr(namespace, self.dest, [v.split(":") for v in values])

def get_arguments():
  parser = ArgumentParser(description="A tool to extract tallies from MCNP output file")
  parser.add_argument(
    '-k', dest='divide_keff',
    help='Divide tallies by keff',
    default=False,
    action='store_true',
  )
  parser.add_argument(
    'files',
    help='File names to process',
    type=str,
    nargs='+',
    action=ColonSplitAction,
    metavar='<file1 [file2 ...]>'
  )
  return parser, parser.parse_args()

def is_float(text):
  try:
    test = float(text)
    return True
  except:
    return False

def is_in_tally(sline, num):
  res = False
  if len(sline) > 2:
    res = True if (sline[0] == '1tally' and sline[1] == num) else False
  return res

def get_tally_numbers(flines):
  tally_list = []
  for line in flines:
    sline = line.split()
    if sline:
      if len(sline) > 2:
        if sline[0] == '1tally' and sline[1] != 'fluctuation':
          tally_list.append(sline[1])
  return tally_list

def get_tally_data(tally_number, flines):
  in_tally, mat, tally_dict, slines = False, None, {}, [l.split() for l in flines]
  for i, line in enumerate(slines):
    in_tally = is_in_tally(line, tally_number) or in_tally
    if in_tally:
      if line:
        if line[0] == '+':
          tally_dict['comment'] = flines[i][1:].strip()
        if line[0] == 'cell' and 'is' not in line:
          cell_number = flines[i][6:].strip()
          if cell_number not in tally_dict.keys():
            tally_dict[cell_number] = {}
          j = 0
          while slines[i+j]:
            if is_float(slines[i+j][0]) or slines[i+j][0] == 'total':
              if len(slines[i+j]) == 2:
                v, fsdev = float(slines[i+j][0]), float(slines[i+j][1])
                if mat:
                  tally_dict[cell_number][mat]['total'] = ufloat(v, v * fsdev)
                else:
                  if '0' not in tally_dict[cell_number]:
                    tally_dict[cell_number]['0'] = {}
                  tally_dict[cell_number]['0']['total'] = ufloat(v, v * fsdev)
              elif len(slines[i+j]) == 3:
                v, fsdev = float(slines[i+j][1]), float(slines[i+j][2])
                if mat:
                  tally_dict[cell_number][mat][slines[i+j][0]] = ufloat(v, v * fsdev)
                else:
                  if '0' not in tally_dict[cell_number]:
                    tally_dict[cell_number]['0'] = {}
                  tally_dict[cell_number]['0'][slines[i+j][0]] = ufloat(v, v * fsdev)
            elif slines[i+j][0] == 'multiplier':
              mat = slines[i+j][3]
              tally_dict[cell_number][mat] = {}
            j += 1
        if line[0] == '*****' or line[0] == 'there':
          return tally_dict
        elif len(line) > 1:
          if line[1] == '=':
            return tally_dict

def get_test_result(tally_number, flines):
  slines = [l.split() for l in flines]
  for i, sline in enumerate(slines):
    if len(sline) == 19:
      if sline[0] == 'results' and sline[12] == 'fluctuation' and sline[-1] == tally_number:
        results = slines[i+7][1:]
        return all(x == 'yes' for x in results)

def get_all_tallies(flines):
  return {tally_num : get_tally_data(tally_num, flines) for tally_num in get_tally_numbers(flines)}

def get_tally_test_results(flines):
  return {tally_num : get_test_result(tally_num, flines) for tally_num in get_tally_numbers(flines)}

if __name__ == '__main__':
  file_dicts, (argp, args) = {}, get_arguments()
  for i, fyle in enumerate(args.files):
    print(f'File: {fyle[0]}')
    keff = get_keff(fyle[0]) if args.divide_keff else 1.
    with open(fyle[0], 'r', errors='ignore') as f:
      lines = f.readlines()
      file_dicts[fyle[0]] = get_all_tallies(lines)
      if len(fyle) > 1:
        new_dict = NestedDefaultDict()
        new_dict[fyle[0]][fyle[1]][fyle[2]][fyle[3]][fyle[4]] = file_dicts[fyle[0]][fyle[1]][fyle[2]][fyle[3]][fyle[4]]
        file_dicts = new_dict
      for tally, cell_data in file_dicts[fyle[0]].items():
        if 'comment' in cell_data:
          c = cell_data['comment']
          print(f'  Tally #: {tally} | {c}')
        else:
          print(f'  Tally #: {tally}')
        for cell, mat_data in cell_data.items():
          if cell != 'comment' and any([any([data.nominal_value > 0. for data in energy_data.values()]) for energy_data in mat_data.values()]):
            print(f'    Cell #: {cell}')
            for mat, energy_data in mat_data.items():
              print(f'      Material #: {mat}')
              for e, data in energy_data.items():
                if data.nominal_value > 0.:
                  if len(args.files) > 1 and i != 0:
                    ref_data = file_dicts[args.files[0]][tally][cell][mat][e]
                    ref_keff = get_keff(args.files[0]) if args.divide_keff else 1.
                    if ref_data.nominal_value > 0.:
                      diff = 100. * ((data/keff) - (ref_data/ref_keff)) / (ref_data/ref_keff)
                      print(f'        energy {e} : {(data/keff):.1uS} [% diff: {diff:.1uS}]')
                  else:
                    print(f'        energy {e} : {(data/keff):.1uS}')
        print(f'  Tests passed for tally #{tally}: {get_test_result(tally, lines)}')
