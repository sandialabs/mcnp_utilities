#! /usr/bin/env python3

from os import walk
from os.path import isdir, isfile, join as path_join
from argparse import ArgumentParser
from uncertainties import ufloat
from pandas import DataFrame


def get_keff(out_file):
  with open(out_file, 'r', errors='ignore') as f:
    for line in f:
      line = line.split()
      if len(line) > 7:
        if line[6] == 'keff':
          return ufloat(float(line[8]), float(line[15]))

def get_beff(out_file):
  with open(out_file, 'r', errors='ignore') as f:
    for line in f:
      line = line.split()
      if len(line) == 3:
        if line[0] == 'beta-eff':
          return ufloat(float(line[1]), float(line[2]))

def get_gen_time(out_file):
  with open(out_file, 'r', errors='ignore') as f:
    for line in f:
      line = line.split()
      if len(line) == 5:
        if line[0] == 'gen.' and line[1] == 'time':
          return ufloat(float(line[2]), float(line[3]))

def get_nubar(out_file):
  with open(out_file, 'r', errors='ignore') as f:
    for line in f:
      line = line.split()
      if len(line) > 11:
        if line[5] == 'neutrons' and line[6] == 'produced' and line[7] == 'per' and line[8] == 'fission':
          return float(line[10])

def calc_rho(keff, beta_eff):
  return (keff - 1.0) / (keff * beta_eff) if keff and beta_eff else None

def get_rho(out_file):
  return calc_rho(get_keff(out_file), get_beff(out_file))

def process_file(out_file, quiet: bool=False) -> dict:
  results = {
    'keff' : get_keff(out_file),
    'beff' : get_beff(out_file)
  }
  results['rho'] = calc_rho(results['keff'], results['beff'])
  if all([x is not None for x in [results['keff'], results['beff'], results['rho']]]):
    if not quiet:
      print(f'{out_file}:')
      for k, v in results.items():
        print("{0:<10} : {1:.1uS}$".format(k, v)) if k == 'rho' else print("{0:<10} : {1:.1uS}".format(k, v))
    return results

def getArguments():
  parser = ArgumentParser(description="A tool to compute the reactivity parameters of an MCNP model",)
  parser.add_argument(
    '-e', dest='extension',
    type=str,
    help='File extension to look for if processing a folder',
    metavar='<extension>',
    default='.out'
  )
  parser.add_argument(
    '-x', dest='excel',
    type=str,
    help='Excel file prefix (without extension) to which to write all results if processing a directory',
    metavar='<excel file>'
  )
  parser.add_argument(
    '-r', dest='reference',
    type=str,
    help='File to compute reactivity relative to',
    metavar='<file>'
  )
  parser.add_argument(
    'path',
    help='Path to output file or folder',
    metavar='<path>'
  )
  return parser.parse_args()

if __name__ == '__main__':
  args = getArguments()
  if args.excel and not isdir(args.path):
    print('Excel file will not be written because path is not a directory!')
  if isfile(args.path):
    res = process_file(args.path)
    if args.reference:
      print(f'Reactivity difference: {(res["rho"]-process_file(args.reference, quiet=True)["rho"]):.1uS}$')
  elif isdir(args.path):
    all_results, ref_rho = {}, process_file(args.reference, quiet=True)['rho'] if args.reference else None
    for root, dirs, files in walk(args.path):
      for _file in files:
        if _file.endswith(args.extension):
          __file = path_join(root, _file)
          res = process_file(__file)
          if res is not None:
            if args.reference:
              print(f'Reactivity difference: {(res["rho"]-ref_rho):.1uS}$')
            all_results[__file] = res
    if args.excel:
      ddict = {k : [] for k in ['File', 'keff', 'beff', 'rho', f'rho diff relative to {args.reference}' if args.reference else None]}
      for i, (f, rdict) in enumerate(all_results.items()):
        for k, v in zip(ddict.keys(), [f, f'{rdict["keff"]:.1uS}', f'{rdict["beff"]:.1uS}', f'{rdict["rho"]:.1uS}', f'{rdict["rho"]-ref_rho:.1uS}' if args.reference else None]):
          ddict[k].append(v)
      DataFrame(ddict).to_excel(f'{args.excel}.xlsx', index=False)
