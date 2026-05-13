#!/usr/bin/env python3

from multiprocessing import Pool, freeze_support
from os import walk, getcwd, remove
from os.path import join as path_join, split as path_split
from argparse import ArgumentParser
# Local modules
from mcnp_utilities.utils.meshtal_to_h5 import write_file


def get_arguments():
  parser = ArgumentParser(description="A tool to convert all meshtal files in a directory (including subdirectories) to hdf5.",)
  parser.add_argument(
    '-f', '--file',
    help='Name of the meshtal files. Default: %(default)s',
    type=str,
    default='meshtal',
    metavar='<file>'
  )
  parser.add_argument(
    '-d','--delete',
    help='Delete the original meshtal files after conversion. Default: %(default)s',
    action='store_true'
  )
  parser.add_argument(
    '-j', '--processes',
    help='Number of threads to use to process files. (default: %(default)s)',
    default=1,
    type=int,
    metavar='<# processes>'
  )
  return parser.parse_args()

def process_file(fpath,delete):
  print("Processing {}...".format(fpath))
  ps = list(path_split(fpath))
  ps[-1] += '.h5'
  out_name = path_join(*ps)
  write_file(out_name,fpath)
  if delete:
    print("Deleting {}...".format(fpath))
    remove(fpath)

def main():
  args = get_arguments()
  mesh_files = [
    path_join(root, name)
    for root, _, files in walk(getcwd())
    for name in files
    if name == args.file
  ]
  with Pool(processes=args.processes) as pool:
    pool.starmap(process_file,[(f,args.delete) for f in mesh_files])

if __name__ == '__main__':
  freeze_support()
  main()
