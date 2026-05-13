#!/usr/bin/env python3

from copy import deepcopy
from os.path import sep, join as pjoin, exists, isdir
from os import mkdir, walk, listdir
from collections import defaultdict
from math import sqrt
import matplotlib.pyplot as plt


def create_nested_path(path: str) -> None:
  spath = proper_path_split(path)
  for l in range(len(spath)):
    p = pjoin(*spath[:l+1])
    if not exists(p):
      mkdir(p)

def write_template_file(file_path, template_lines, data_dict, start='{', end='}'):
  with open(file_path, 'w') as f:
    for line in template_lines:
      for k, v in data_dict.items():
        line = line.replace(f'{start}{k}{end}', str(v)) if f'{start}{k}{end}' in line else line
      f.write(line)

def set_plot_params():
  params = {
    "text.usetex" : True,
    "font.family" : "serif",
    "font.serif" : ["Computer Modern Serif"]
  }
  plt.rcParams.update(params)

def proper_path_split(path: str) -> list[str]:
  return [x for x in path.split(sep) if x]

def has_subdirectory(path: str) -> bool:
  return any(isdir(pjoin(path, x)) for x in listdir(path))

class bottomlessdict(defaultdict):
  def __init__(self) -> None:
    self.default_factory = lambda: defaultdict(defaultdict)

def get_directorystructuredict(path: str, exclusions=[], terminator=None) -> dict:
  this_d = {}
  for root, directories, _ in walk(path, topdown=True):
    for directory in directories:
      d, keys, cwdkeys = this_d, proper_path_split(pjoin(root, directory)), proper_path_split(path)
      for idx, folder in enumerate(keys):
        if (idx >= len(cwdkeys)) and (folder not in exclusions):
          if folder not in d:
            d[folder] = dict() if has_subdirectory(pjoin(root, directory)) else deepcopy(terminator)
          d = d[folder]
  return this_d

def hex_2_cart(x: float, y: float, P: float):
  x_cart = P * (x + y / 2)
  y_cart = P * y * sqrt(3) / 2
  return x_cart, y_cart

def write_modified_file(original_file: str, data_dict: dict, new_file: str=None, start='{', end='}') -> None:
  with open(original_file, 'r', errors='ignore') as f:
    original_lines = f.readlines()
  destination_file = new_file if new_file is not None else original_file
  write_template_file(destination_file, original_lines, data_dict, start=start, end=end)
