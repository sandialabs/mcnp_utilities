#!/usr/bin/env python3

import os

files = [f for f in os.listdir('.') if os.path.isfile(f)]
rm_files = ['out', 'runtp', 'srct', 'comou']
for file in files:
  if any(file[0:len(x)] == x for x in rm_files):
    os.remove(file)
