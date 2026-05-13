#!/usr/bin/env python3

"""
Steadystate is a tool that allows for manual and automated searching
for target eigenvalue cases using MCNP by changing the value of a
variable defined in the MCNP input file.

The manual search is the same as running several different MCNP files
and comparing the results. On the other hand, using the automated
bisection search will allow for more precision and less computational
work.

# Author    : Zackary Dodson
"""

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from os import getcwd, makedirs, devnull, remove, listdir, rename, walk
from os.path import exists, isfile, join
from shutil import move, rmtree
from subprocess import check_call
from uncertainties import ufloat
from math import log, exp


# =============================================================================
# Globals
pwd = getcwd()
steadystate_title_text = r"""      ____________________________________________________
     / ______              __     ______       __        /
    / / __/ /____ ___ ____/ /_ __/ __/ /____ _/ /____   /
   / _\ \/ __/ -_) _ `/ _  / // /\ \/ __/ _ `/ __/ -_) /
  / /___/\__/\__/\_,_/\_,_/\_, /___/\__/\_,_/\__/\__/ /
 /                        /___/                      /
/___________________________________________________/
  """
# -----------------------------------------------------------------------------
def Steadystate_Input_Parser():
  parser = ArgumentParser(
    description='This program finds the closest value for a given keff or reactivity target.',
    formatter_class=ArgumentDefaultsHelpFormatter
  )
  run_type = parser.add_mutually_exclusive_group()
  run_type.add_argument(
    '-m', '--manual',
    help='File containing newline-delimited variable values to use for search.',
    metavar='<file name>'
  )
  run_type.add_argument(
    '-r', '--ranged',
    nargs=3,
    type=float,
    metavar=('<lowest value>', '<range>', '<# of intervals>')
  )
  run_type.add_argument(
    '-b', '--bisection',
    nargs=4,
    type=float,
    metavar=('<bound 1>', '<bound 2>', '<max iterations>', '<tolerance (pcm|$)>')
  )
  run_type.add_argument(
    '-l', '--linear-regression',
    nargs=4,
    type=float,
    metavar=('<bound 1>', '<bound 2>', '<desired unc. (pcm|$)>', '<tolerance (pcm|$)>')
  )
  parser.add_argument(
    '-t', '--tasks',
    help='Number of processors to use for MCNP calculations. Default: %(default)s',
    default=1,
    metavar='<# of threads>'
  )
  parser.add_argument(
    '-s', '--source',
    help='Source tape file to use for MCNP calculations.',
    metavar='<file name>'
  )
  parser.add_argument(
    '-d', '--decimals',
    help='Number of decimal places to round variable value before it is inserted into the MCNP input file. Default: %(default)s',
    type=int,
    default=15,
    metavar='<# of decimals>'
  )
  parser.add_argument(
    '-e', '--executable',
    help='Path to MCNP executable. If no value is entered, it is assumed that the path is defined in the PATH environmental variable. Default: %(default)s',
    default='mcnp6',
    metavar='<path>'
  )
  parser.add_argument(
    '-impi', dest='impi',
    help='Used when using Intel MPI. Argument is # of processes per node.',
    type=int,
    metavar='<# ppn>'
  )
  parser.add_argument(
    '--np',
    help='Number of processors to use for mpi calculation.',
    metavar='<np>'
  )
  parser.add_argument(
    '-k', '--keep',
    help='Flag for keeping input and output files.',
    action='store_true'
  )
  parser.add_argument(
    '-c', '--clean',
    help='Flag for clearing all files and directories made by this program in the current directory.',
    action='store_true'
  )
  parser.add_argument(
    '-n', '--mcnpclean',
    help='Flag for deleting all MCNP-related files (starting with runtp, '
         'srct, out, comou, and meshta as well as MCNP input files generated '
         'by this script). WARNING: This includes any file containing '
         '"_st_" in the name.',
    action='store_true'
  )
  parser.add_argument(
    '-q', '--quiet',
    help='Suppresses MCNP output. Makes it easy to keep track of progress.',
    action='store_true'
  )
  parser.add_argument(
    'file',
    type=str,
    help='File to process',
    metavar='<file name>'
  )
  parser.add_argument(
    'target',
    type=float,
    help='Target keff or reactivity (in dollars) for bisection or linear regression case',
    metavar='<target>'
  )
  parser.add_argument(
    'type',
    type=str,
    help='Specifies the type of target (keff or rho)',
    choices=['keff', 'rho']
  )
  return parser.parse_args()
# -----------------------------------------------------------------------------
def Steadystate(args):
  sign = lambda x: (1, -1)[x < 0]
  if args.clean:
    clean()
  if args.mcnpclean:
    MCNP_Clean()
  if not args.ranged and not args.bisection and not args.manual and not args.linear_regression:
    raise RuntimeError('ERROR: Input is incorrect. Missing ranged, manual, bisection, or linear regression arguments.')
  input_path = join(pwd, 'Input_Files')
  output_path = join(pwd, 'Output_Files')
  if args.keep:
    if not exists(input_path):
      makedirs(input_path)
    if not exists(output_path):
      makedirs(output_path)
  if not args.quiet:
    print('File to be run: %s\n' % args.file)
  loc = args.file.find('.')
  base_file_name = args.file[:loc] if loc > -1 else args.file
  paths = [pwd, input_path, output_path]
  lsts = {s : [] for s in ['values', 'error', 'pk_vals']}
  vloclist, mem = get_data(args.file)
  case = 0

  # Handle ranged cases
  if args.ranged or args.manual:
    if args.ranged:
      step, cases = get_ranged_cases(args.ranged[0], args.ranged[1], args.ranged[2], args.decimals)
    else:
      cases = get_manual_cases(args.manual)

    # Loop over all cases
    for value in cases:
      pk_vals, abserr, err, case, problem = run_case(base_file_name, case, value, vloclist, mem, args, lsts, paths)
      if problem:
        break

  # Handle bisection cases
  if args.bisection:
    case_limit, err_tol = int(args.bisection[2]), args.bisection[3]
    min_val = min([args.bisection[0], args.bisection[1]])
    max_val = max([args.bisection[0], args.bisection[1]])
    endpts = [min_val, max_val]
    if not args.quiet:
      print('Performing search using bisection method.')
      if args.type == 'rho':
        print('Target rx: ${}'.format(args.target))
      elif args.type == 'keff':
        print('Target k-eff: {}'.format(args.target))
      print('Lower bound : {}'.format(min_val))
      print('Upper bound : {}'.format(max_val))
      print('A maximum of {0:.0f} iterations will be allowed.'.format(case_limit))
      if args.type == 'keff':
        print('A tolerance of {0:.0f} pcm is allowed.'.format(err_tol))
      elif args.type == 'rho':
        print('A tolerance of ${0:.3f} is allowed.'.format(err_tol))

    endpt_abserr, problem = [0, 0], False
    while True:
      if case > case_limit:
        print('The specified case limit of {0} has been reached without satisfying the termination criterion.'.format(case_limit))
        break
      if case <= 1:
        for i, endpt in enumerate(endpts):
          pk_vals, endpt_abserr[i], err, case, problem = run_case(base_file_name, case, endpt, vloclist, mem, args, lsts, paths, err_tol)
          if abs(err) < err_tol:
            if args.type == 'keff':
              print('Endpoint {0} is within the specified error tolerance of {1:.0f} pcm with a pcm error of {2} and keff = {3}.\n'.format(endpt, err_tol, err, pk_vals['k'].nominal_value))
            elif args.type == 'rho':
              print('Endpoint {0} is within the specified error tolerance of ${1:.2f} with a error of ${2} and rho = ${3}.\n'.format(endpt, err_tol, err, pk_vals['rho'].nominal_value))
            done = True
            break
        else:
          if sign(endpt_abserr[0]) == sign(endpt_abserr[1]):
            raise RuntimeError('The target value does not appear to lie between the specified points.')
      else:
        binary_iterate(min_val, max_val, args, base_file_name, case, vloclist, mem, lsts, paths, err_tol, case_limit)
        break

  if args.linear_regression:
    run_linear_regression(args, base_file_name, vloclist, mem, paths, lsts)

  # Find the closest case to target.
  if args.type == 'keff':
    pos, min_index = get_closest_position(lsts["values"], [lsts["pk_vals"][i]['k'] for i in range(len(lsts["pk_vals"]))], args.target)
    print('The closest value for {0} to the target k-eff of {1:.5f} is: {2}.'.format(args.file, args.target, pos))
  elif args.type == 'rho':
    pos, min_index = get_closest_position(lsts["values"], [lsts["pk_vals"][i]['rho'] for i in range(len(lsts["pk_vals"]))], args.target)
    print('The closest value for {0} to the target rho of ${1:.5f} is: {2}.'.format(args.file, args.target, pos))
  if __name__ == '__main__':
    # Print results.
    titles = []
    text = make_results(lsts, titles, args.target, args.type)
    # Write results to output.
    write_results(base_file_name, lsts, titles, text, args.target)
    print('''\nAll runs complete.''')
  else:
    return pos, lsts["pk_vals"][min_index]
# -----------------------------------------------------------------------------
def run_linear_regression(args, base_file_name, vloclist, mem, paths, lsts):
  """
  This method is an implementation of the GRsecant method in the following reference:

  "Method for Control Drum Position Critical Search with Monte Carlo Codes" by Dean Price & Nathan Roskoff
  """
  # Run endpoints and compute distances to target
  pk_vals, abserr, err, case, problem = run_case(base_file_name, 0, args.linear_regression[0], vloclist, mem, args, lsts, paths, args.linear_regression[3])
  pk_vals, abserr, err, case, problem = run_case(base_file_name, 1, args.linear_regression[1], vloclist, mem, args, lsts, paths, args.linear_regression[3])
  if args.type == 'keff':
    f = [(lsts['pk_vals'][0]['k'] - args.target).nominal_value, (lsts['pk_vals'][1]['k'] - args.target).nominal_value]
  elif args.type == 'rho':
    f = [(lsts['pk_vals'][0]['rho'] - args.target).nominal_value, (lsts['pk_vals'][1]['rho'] - args.target).nominal_value]
  # Setup
  n, R = 2, 2
  g0 = get_initial_active_cycles(mem)
  g = [g0, g0]
  # Iterate
  while True:
    # Compute new value and round
    val = round(compute_linear_regression_value(lsts, R, f, args), args.decimals)
    # Ensure value is in bounds
    val = max(args.linear_regression[0], val)
    val = min(args.linear_regression[1], val)
    # Compute target uncertainty
    if args.type == 'keff':
      tgt_unc = 0.95 * (args.linear_regression[2] / 1e5) * (min([abs(x) for x in f]) / (args.linear_regression[3] / 1e5))**(0.5)
      if tgt_unc < 0.95 * (args.linear_regression[2] / 1e5):
        tgt_unc = 0.95 * (args.linear_regression[2] / 1e5)
    elif args.type == 'rho':
      tgt_unc = 0.95 * args.linear_regression[2] * (min([abs(x) for x in f]) / args.linear_regression[3])**(0.5)
      if tgt_unc < 0.95 * args.linear_regression[2]:
        tgt_unc = 0.95 * args.linear_regression[2]
    # Compute new number of active generations to run and update input
    compute_new_active_gens(tgt_unc, lsts, g, n, args)
    update_active_gens(mem, g[-1])
    # Run case and compute distance from target
    pk_vals, abserr, err, n, problem = run_case(base_file_name, n, val, vloclist, mem, args, lsts, paths, args.linear_regression[3])
    if args.type == 'keff':
      f.append((lsts['pk_vals'][-1]['k'] - args.target).nominal_value)
    elif args.type == 'rho':
      f.append((lsts['pk_vals'][-1]['rho'] - args.target).nominal_value)
    # Check for convergence
    if args.type == 'keff':
      if abs(err) < args.linear_regression[3] and lsts['pk_vals'][-1]['k'].std_dev < (args.linear_regression[2] / 1e5):
        break
    elif args.type == 'rho':
      if abs(err) < args.linear_regression[3] and lsts['pk_vals'][-1]['rho'].std_dev < args.linear_regression[2]:
        break
# -----------------------------------------------------------------------------
def compute_linear_regression_value(lsts, R, f, args):
  s1, s2, s3, s4, s5 = 0, 0, 0, 0, 0
  for i in range(R):
    if args.type == 'keff':
      v = lsts["pk_vals"][-1 - i]['k']
    elif args.type == 'rho':
      v = lsts["pk_vals"][-1 - i]['rho']
    c = lsts['values'][-1 - i] / (v.std_dev)**2
    s1 += c
    s2 += f[-1 - i] * c
    s3 += (lsts['values'][-1 - i])**2 / (v.std_dev)**2
    s4 += f[-1 - i] / (v.std_dev)**2
    s5 += 1.0 / (v.std_dev)**2
  return (s1 * s2 - s3 * s4) / (s2 * s5 - s4 * s1)
# -----------------------------------------------------------------------------
def compute_new_active_gens(tgt_unc, lsts, g, n, args):
  s1, s2, s3, s4 = 0, 0, 0, 0
  for i in range(n):
    if args.type == 'keff':
      v = lsts["pk_vals"][-1 - i]['k']
    elif args.type == 'rho':
      v = lsts["pk_vals"][-1 - i]['rho']
    s1 += log(v.std_dev)
    s2 += log((g[i]))**2
    s3 += log(v.std_dev) * log(g[i])
    s4 += log(g[i])
  t1 = ((s1 * s2 - s3 * s4) / ((n + 1) * s2 - ((s4)**2)))
  t2 = (((n + 1) * s2 - ((s4)**2)) / (s1 * s4 - (n + 1) * s3))
  g.append(max(30, round(exp((t1 - log(tgt_unc)) * t2))))
# -----------------------------------------------------------------------------
def binary_iterate(min_val, max_val, args, base_file_name, case, vloclist, mem, lsts, paths, err_tol, case_limit):
  mid = round((min_val + max_val) / 2, args.decimals)
  pk_vals, mid_err, err, case, problem = run_case(base_file_name, case, mid, vloclist, mem, args, lsts, paths, err_tol)
  v = pk_vals['k'].nominal_value if args.type == 'keff' else pk_vals['rho'].nominal_value
  if abs(err) < err_tol or case > case_limit:
    return mid
  elif v < args.target:
    return binary_iterate(mid, max_val, args, base_file_name, case, vloclist, mem, lsts, paths, err_tol, case_limit)
  else:
    return binary_iterate(min_val, mid, args, base_file_name, case, vloclist, mem, lsts, paths, err_tol, case_limit)
# -----------------------------------------------------------------------------
def write_results(fname, lsts, titles, text, target):
  with open(fname + '_results.txt', 'w') as f:
    [f.write(line + '\n') for line in text]
# -----------------------------------------------------------------------------
def make_results(lsts, titles, target, ttype, quiet=False):
  text = []
  has_beta = True if None not in [lsts["pk_vals"][i]['b'].nominal_value for i in range(len(lsts["pk_vals"]))] else False
  if ttype == 'keff':
    titles += [' Value', ' keff', ' keff-sdev', ' pcm-sdev', ' pcm-error']
    title_line = '{0:10}|{1:10}|{2:10}|{3:10}|{4:10}'.format(titles[0], titles[1], titles[2], titles[3], titles[4])
  elif ttype == 'rho':
    titles += [' Value', ' rho ', ' rho-sdev ', ' rho-error']
    title_line = '{0:10}|{1:10}|{2:10}|{3:10}'.format(titles[0], titles[1], titles[2], titles[3])
  if has_beta:
    titles += [' Beff', ' Beff-sdev']
    if ttype == 'keff':
      title_line += '|{0:10}|{1:10}'.format(titles[5], titles[6])
    elif ttype == 'rho':
      title_line += '|{0:10}|{1:10}'.format(titles[4], titles[5])
  text.append(title_line)
  if not quiet:
    print(title_line)
  for i in range(len(lsts["values"])):
    if ttype == 'keff':
      pcm_sdev = int(round(lsts["pk_vals"][i]['k'].std_dev * 1E+05, 0))
      row = '{0:< 10.8}|{1:< 10.6}|{2:< 10.6}|{3:< 10}|{4:< 10}'.format(float(lsts["values"][i]), lsts["pk_vals"][i]['k'].nominal_value, lsts["pk_vals"][i]['k'].std_dev, pcm_sdev, lsts["error"][i])
    elif ttype == 'rho':
      row = '{0:< 10.8}|{1:< 10.6}|{2:< 10.6}|{3:< 10}'.format(float(lsts["values"][i]), lsts["pk_vals"][i]['rho'].nominal_value, lsts["pk_vals"][i]['rho'].std_dev, lsts["error"][i])
    if has_beta:
      row += '|{0:< 10}|{1:< 10}'.format(lsts["pk_vals"][i]['b'].nominal_value, lsts["pk_vals"][i]['b'].std_dev)
    if not quiet:
      print(row)
    text.append(row)
  return text
# -----------------------------------------------------------------------------
def get_closest_position(value_list, k_list, target):
  if len(value_list) != len(k_list):
    raise AssertionError('List lengths are not the same!')
  chklist = [abs(target - k) for k in k_list]
  min_index = chklist.index(min(chklist))
  return value_list[min_index], min_index
# -----------------------------------------------------------------------------
def get_manual_cases(manual_inp):
  if isinstance(manual_inp, str):
    cases = []
    with open(manual_inp, 'r') as f:
      [cases.append(float(line)) for line in f]
  else:
    cases = manual_inp
  return cases
# -----------------------------------------------------------------------------
def get_ranged_cases(low_val, rng, ints, rnd, quiet=None):
  case, cases = 0, []
  step = rng / ints if ints > 0 else 0
  # Create case list
  if not quiet:
    print('Variable values to be used:')
  while case <= ints:
    val = round(low_val + (step * case), rnd)
    cases.append(val)
    if not quiet:
      print(str(val))
    case += 1
  return step, cases
# -----------------------------------------------------------------------------
def get_data(input_file):
  """
  Finds the variable string designator ({X}) in the input file and
  records the location. Also stores all of the input file into a list.

  input_file : The MCNP input file to be run.
  """
  i, found_flag, mem, vloclist = 0, False, [], []
  with open(str(input_file), 'r') as f:
    for i, line in enumerate(f):
      mem.append(line)
      try:
        if '{X}' in line:
          vloclist.append(i)
          found_flag = True
      except:
        pass
    if not found_flag:
      print('''\nError: {X} not found.''')
      quit()
  return vloclist, mem
# -----------------------------------------------------------------------------
def variable_change(newval, vloclist, mem, name, args):
  """
  Changes the value of the variable in the input file.

  newval   : The new value of the variable.
  vloclist : List of variable line number locations.
  mem      : List of lines in input file.
  name     : New file name (with {X} replaced with number).
  args     : Input argument list.
  """
  with open(name, 'w') as f:
    for i, line in enumerate(mem):
      write_line = line
      if i in vloclist:
        write_line = write_line.replace('{X}', str(round(newval, args.decimals)))
      f.write(write_line)
# -----------------------------------------------------------------------------
def get_initial_active_cycles(mem):
  for i, line in enumerate(mem):
    sline = line.lower().split()
    if len(sline) > 4:
      if sline[0] == 'kcode':
        return int(sline[4]) - int(sline[3])
# -----------------------------------------------------------------------------
def update_active_gens(mem, gens):
  for i, line in enumerate(mem):
    sline = line.lower().split()
    if len(sline) > 4:
      if sline[0] == 'kcode':
        sline[4] = str(max(30, gens) + int(sline[3]))
        mem[i] = " ".join(sline) + '\n'
        break
# -----------------------------------------------------------------------------
def run_MCNP(input_file_name, case, base_file_name, paths, args):
  """
  Runs MCNP with the current input file state, handles the output
  files, and finds keff in the output.

  input_file_name : Name of the input file.
  case            : Case number.
  base_file_name  : Base file name (no extension or case number
                    specifier).
  paths           : File paths.
  args            : Input argument list.
  """
  out_name = base_file_name + '_case_' + str(case) + '.out'
  if not args.np:
    if args.source:
      com = ('%s i=%s o=%s s=%s tasks %s\n' % (args.executable, input_file_name, out_name, args.source, args.tasks))
    else:
      com = ('%s i=%s o=%s tasks %s\n' % (args.executable, input_file_name, out_name, args.tasks))
  else:
    if args.source:
      if args.impi:
        com = ('mpiexec -ppn %s --bind-to socket -n %s %s i=%s o=%s s=%s tasks %s\n' % (args.impi, args.np, args.executable, input_file_name, out_name, args.source, args.tasks))
      else:
        com = ('mpiexec --bind-to socket --npersocket 1 --n %s %s i=%s o=%s s=%s tasks %s\n' % (args.np, args.executable, input_file_name, out_name, args.source, args.tasks))
    else:
      if args.impi:
        com = ('mpiexec -ppn %s --bind-to socket -n %s %s i=%s o=%s tasks %s\n' % (args.impi, args.np, args.executable, input_file_name, out_name, args.tasks))
      else:
        com = ('mpiexec --bind-to socket --npersocket 1 --n %s %s i=%s o=%s tasks %s\n' % (args.np, args.executable, input_file_name, out_name, args.tasks))
  if args.quiet:
    FNULL = open(devnull, 'w')
    check_call(com, shell=True, stdout=FNULL)
  else:
    check_call(com, shell=True)
  if isfile('runtpe.h5'):
    remove('runtpe.h5')
  if isfile('runtpe'):
    remove('runtpe')
  if isfile('srctp'):
    remove('srctp')
  if args.source:
    [remove(f) for f in listdir('.') if isfile(f) and
     f.startswith(args.source[:-1]) and f != args.source]
  if isfile('meshtal'):
    meshname = ''.join([base_file_name,'_case_',str(case),'.mesh'])
    rename('meshtal', meshname)
    meshpath = join(paths[0], 'Mesh')
    if not exists(meshpath):
      makedirs(meshpath)
    move(join(paths[0], meshname, meshpath, meshname))
  pk_vals, problem = parse_output(out_name)
  if problem:
    pk_vals['k'].nominal_value, pk_vals['b'].std_dev = -1.0, -1.0
  if not args.keep:
    remove(out_name)
  if isfile(input_file_name) and not args.keep:
    remove(input_file_name)
  else:
    move(join(paths[0], input_file_name), join(paths[1], input_file_name))
    move(join(paths[0], out_name), join(paths[2], out_name))
  return pk_vals, out_name, problem
# -----------------------------------------------------------------------------
def parse_output(file):
  """
  Obtain point kinetics parameters from MCNP output file.

  file : file to be parsed
  """
  problem = True
  pk_vals = {}
  with open(file, 'r', errors='ignore') as f:
    for line in f:
      line = line.split()
      if len(line) > 6 and line[6] == 'keff':
        problem = False
        pk_vals['k'] = ufloat(float(line[8]), float(line[15]))
      elif len(line) == 3 and line[0] == 'beta-eff':
        pk_vals['b'] = ufloat(float(line[1]), float(line[2]))
      elif len(line) == 5 and line[0] == 'gen.' and line[1] == 'time':
        pk_vals['l'] = ufloat(float(line[2]), float(line[3]))
  if pk_vals['k'] and pk_vals['b']:
    pk_vals['rho'] = (pk_vals['k'] - 1.0) / (pk_vals['k'] * pk_vals['b'])
  return pk_vals, problem
# -----------------------------------------------------------------------------
def run_case(base_file_name, case, newval, vloclist, mem, args, lsts, paths, e_tol=None):
  """
  Handle a single case. Calls functions to change variable value and
  run the new file in MCNP. Then extracts output data.

  base_file_name  : Base file name (no extension or case number specifier).
  case            : Case number.
  newval          : New variable value.
  vloclist        : List of variable line number locations.
  mem             : List of lines in input file.
  args            : Input argument list.
  lsts            : Case data dictionary.
  paths           : File paths.
  """
  input_file_name = ''.join([base_file_name,'_case_',str(case),'.inp'])
  variable_change(newval, vloclist, mem, input_file_name, args)
  lsts["values"].append(newval)
  print('~' * 79)
  print('Running case {0} for file {1} with value {2}.'.format(case, input_file_name, newval))
  pk_vals, out_name, problem = run_MCNP(input_file_name, case, base_file_name, paths, args)
  if args.type == 'keff':
    err = int(round((pk_vals['k'].nominal_value - args.target) * 100000, 0))
    print('keff = {0:.5f} | sdev = {1:.5f} | pcm error = {2} | value = {3} '.format(pk_vals['k'].nominal_value, pk_vals['k'].std_dev, err, newval))
    abserr = pk_vals['k'].nominal_value - args.target
  elif args.type == 'rho':
    err = round((pk_vals['rho'].nominal_value - args.target), 3)
    print('rho = ${0:.3f} | sdev = ${1:.3f} | rx error = ${2} | value = {3} '.format(pk_vals['rho'].nominal_value, pk_vals['rho'].std_dev, err, newval))
    abserr = pk_vals['rho'].nominal_value - args.target
  lsts["error"].append(err)
  print('~' * 79)
  lsts["pk_vals"].append(pk_vals)
  return pk_vals, abserr, err, case + 1, problem
# -----------------------------------------------------------------------------
def clean():
  """
  Deletes folders (and contents) created by this script.
  """
  for p in ['Input_Files', 'Output_Files', 'Results']:
    if exists(join(pwd, p)):
      rmtree(join(pwd, p))
# -----------------------------------------------------------------------------
def MCNP_Clean():
  """
  Deletes all MCNP output files and any MCNP input files generated by
  this program.
  """
  for file in next(walk(getcwd()))[2]:
    if any([file.startswith(prefix) for prefix in ['comou', 'out', 'runtp', 'srct', 'meshta']]) or any([s in file for s in ['_st_', '_case_']]):
      remove(file)
# -----------------------------------------------------------------------------
if __name__ == '__main__':
  print(steadystate_title_text)
  Steadystate(Steadystate_Input_Parser())
# -----------------------------------------------------------------------------
