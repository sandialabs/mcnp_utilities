#!/usr/bin/env python3

from os import linesep
from os.path import dirname, join as pjoin
from yaml import safe_load
from collections import defaultdict
# Local modules
from mcnp_utilities.lib.tables import PrintedTable, to_latex_scientific_notation


with open(pjoin(dirname(__file__), '..', 'data', 'iso_ref.yaml'), 'r') as f:
  iso_ref = safe_load(f)

def zaid_to_ZA(zaid_str) -> tuple[int, int]:
  return int(zaid_str[:-3]), int(zaid_str[-3:])

def ZA_to_zaid(Z, A) -> str:
  return f'{Z}{str(A).zfill(3)}'

def get_elemental_zaid(Z) -> str:
  return ZA_to_zaid(Z, 0)

class Isotope:
  def __init__(self, id_str):
    if id_str.isnumeric():
      self.Z, self.A = zaid_to_ZA(id_str)
      if self.Z in iso_ref:
        if self.A not in iso_ref[self.Z]['isotopes']:
          raise AssertionError(f'Isotope {self.Z}-{self.A} not in database!')
        else:
          self.name = f'{iso_ref[self.Z]["name"]}-{self.A}'
      else:
        raise AssertionError(f'Element {self.Z} not in database!')
    else:
      symbol_or_name, self.A = id_str.split('-')[0], int(id_str.split('-')[1])
      for z in iso_ref:
        if iso_ref[z]['name'].lower() == symbol_or_name.lower() or iso_ref[z]['symbol'].lower() == symbol_or_name.lower():
          if self.A in iso_ref[z]['isotopes']:
            self.Z, self.name = int(z), f'{iso_ref[z]["name"]}-{self.A}'
            break
          else:
            raise AssertionError(f'Isotope {symbol_or_name}-{self.A} not in database!')
      else:
        raise AssertionError(f'Element {symbol_or_name} not in database!')
    self.zaid = ZA_to_zaid(self.Z, self.A)

  def get_mass(self):
    Z, A = zaid_to_ZA(self.zaid)
    return iso_ref[Z]['isotopes'][A]['atomic mass']

  def __repr__(self):
    return self.name

class Element:
  def __init__(self, Z, atom_fractions: dict[int, float]={}, weight_fractions: dict[int, float]={}, exclusion_list: list[int]=[]):
    self.Z: int = Z
    self.atom_fractions, self.weight_fractions, self.isotope_list = {}, {}, []

    if atom_fractions and weight_fractions:
      raise AssertionError('Cannot pass both atom and weight fractions to Element initialization!')

    if self.Z not in iso_ref:
      raise AssertionError(f'Element {self.Z} not in database!')
    else:
      self.name = iso_ref[self.Z]['name']

    self.zaid = ZA_to_zaid(self.Z, 0)
    self.exclusion_list = exclusion_list

    if exclusion_list:
      self.exclusion_list = exclusion_list
      self.name = f'{self.name} (no {",".join(str(a) for a in self.exclusion_list)})'

    if atom_fractions:
      self.atom_fractions = atom_fractions
      self.isotope_list = list(self.atom_fractions.keys())
    elif weight_fractions:
      self.weight_fractions = weight_fractions
      self.isotope_list = list(self.weight_fractions.keys())
    else:
      for A in get_natural_isotopes(Z):
        if A not in exclusion_list:
          self.atom_fractions[A] = get_atom_fraction(ZA_to_zaid(self.Z, A))
          self.weight_fractions[A] = get_mass_fraction(ZA_to_zaid(self.Z, A))
          self.isotope_list.append(A)

    self.normalize_fractions()

  def normalize_fractions(self):
    wgt_dict = {ZA_to_zaid(self.Z, A) : get_isotope_mass(ZA_to_zaid(self.Z, A)) for A in self.isotope_list}

    if not self.weight_fractions:
      afs = {ZA_to_zaid(self.Z, k) : v for k, v in self.atom_fractions.items()}
      self.weight_fractions = {k : af_to_wf(ZA_to_zaid(self.Z, k), afs, wgt_dict) for k in self.atom_fractions}
    elif not self.atom_fractions:
      wfs = {ZA_to_zaid(self.Z, k) : v for k, v in self.weight_fractions.items()}
      self.atom_fractions = {k : wf_to_af(ZA_to_zaid(self.Z, k), wfs, wgt_dict) for k in self.weight_fractions}

    af_sum = sum(self.atom_fractions.values())
    wf_sum = sum(self.weight_fractions.values())

    self.atom_fractions = {k : (v / af_sum) for k, v in self.atom_fractions.items()}
    self.weight_fractions = {k : (v / wf_sum) for k, v in self.weight_fractions.items()}

  def get_mass(self):
    return sum([af * get_isotope_mass(ZA_to_zaid(self.Z, A)) for A, af in self.atom_fractions.items()])

  def __repr__(self):
    return self.name

def get_natural_isotopes(Z):
  return [k for k in iso_ref[Z]['isotopes'] if 'atom fraction' in iso_ref[Z]['isotopes'][k]]

def get_mass_fraction(zaid):
  Z, A = zaid_to_ZA(zaid)
  return iso_ref[Z]['isotopes'][A]['mass fraction']

def get_atom_fraction(zaid):
  Z, A = zaid_to_ZA(zaid)
  return iso_ref[Z]['isotopes'][A]['atom fraction']

def get_isotope_mass(zaid):
  Z, A = zaid_to_ZA(zaid)
  return iso_ref[Z]['isotopes'][A]['atomic mass']

def wf_to_af(iso, wf_dict, weight_dict):
  return (wf_dict[iso] / weight_dict[iso]) / (sum([wf_dict[isotope] / weight_dict[isotope] for isotope in wf_dict]))

def af_to_wf(iso, af_dict, weight_dict):
  return (af_dict[iso] * weight_dict[iso]) / (sum([af_dict[isotope] * weight_dict[isotope] for isotope in af_dict]))

def str_to_Z(symbol_or_name):
  for z in iso_ref:
    if iso_ref[z]['name'].lower() == symbol_or_name.lower() or iso_ref[z]['symbol'].lower() == symbol_or_name.lower():
      return z
  else:
    raise AssertionError(f'Element {symbol_or_name} not in database!')

def nuclide_str_to_nuclide(nuc_str):
  if '-' in nuc_str:
    s, A = nuc_str.split('-')
    if not s.isnumeric():
      Z = str_to_Z(s)
    else:
      Z = int(s)
    if int(A) != 0:
      return Isotope(ZA_to_zaid(Z, A))
    else:
      return Element(str_to_Z(nuc_str.split('{')[0]), exclusion_list=[int(a) for a in nuc_str.split('{')[1].split('}')[0].split('+')]) if '{' and '}' in nuc_str else Element(Z)
  else:
    if nuc_str.isnumeric():
      if nuc_str.endswith('000'):
        return Element(zaid_to_ZA(nuc_str.split('{')[0])[0], exclusion_list=[int(a) for a in nuc_str.split('{')[1].split('}')[0].split('+')]) if '{' and '}' in nuc_str else Element(zaid_to_ZA(nuc_str)[0])
      else:
        return Isotope(nuc_str)
    else:
      return Element(str_to_Z(nuc_str.split('{')[0]), exclusion_list=[int(a) for a in nuc_str.split('{')[1].split('}')[0].split('+')]) if '{' and '}' in nuc_str else Element(str_to_Z(nuc_str))

class Material:
  def __init__(self, frac_type, nuclides=[], number=1, density=None, name=None):
    self.number = number
    self.density = density
    self.name = name
    self.nuclides = []
    self.nuclide_atom_fractions, self.nuclide_weight_fractions = defaultdict(float), defaultdict(float)
    if frac_type.lower() == 'atom' or frac_type.lower() == 'weight':
      self.frac_type = frac_type.lower()
    else:
      raise ValueError(f'Unrecognized fraction type "{frac_type}"!')
    if nuclides:
      for n, f in nuclides.items():
        self.add_nuclide(nuclide_str_to_nuclide(n), f)
    self.normalize()

  def add_nuclide(self, nuclide, fraction):
    # Append nuclide to nuclide list if not already present
    if nuclide.zaid not in [n.zaid for n in self.nuclides]:
      self.nuclides.append(nuclide)

    # Add fraction to the correct dictionary
    if self.frac_type == 'atom':
      self.nuclide_atom_fractions[nuclide.zaid] += fraction
    elif self.frac_type == 'weight':
      self.nuclide_weight_fractions[nuclide.zaid] += fraction

  def split_elements(self, split_carbon=True):
    # Create new list of nuclides
    new_nuclide_list = []

    # Loop over nuclides and sum atom or weight fractions while appending to new nuclide list and deleting element entries from fraction dictionaries
    for nuclide in self.nuclides:
      if isinstance(nuclide, Element) and (nuclide.Z != 6 or split_carbon):
        for A in nuclide.isotope_list:
          zaid = ZA_to_zaid(nuclide.Z, A)
          if self.frac_type == 'atom':
            if zaid not in self.nuclide_atom_fractions:
              new_nuclide_list.append(Isotope(zaid))
            self.nuclide_atom_fractions[zaid] += self.nuclide_atom_fractions[nuclide.zaid] * nuclide.atom_fractions[A]
          elif self.frac_type == 'weight':
            if zaid not in self.nuclide_weight_fractions:
              new_nuclide_list.append(Isotope(zaid))
            self.nuclide_weight_fractions[zaid] += self.nuclide_weight_fractions[nuclide.zaid] * nuclide.weight_fractions[A]
        if nuclide.zaid in self.nuclide_atom_fractions:
          del self.nuclide_atom_fractions[nuclide.zaid]
        if nuclide.zaid in self.nuclide_weight_fractions:
          del self.nuclide_weight_fractions[nuclide.zaid]
      else:
        new_nuclide_list.append(nuclide)

    # Set new nuclide list and normalize isotope fractions
    self.nuclides = new_nuclide_list
    self.normalize()

  def combine_isotopes(self):
    new_nuclide_list, fracs, rm_zaids = [], defaultdict(float), []

    # Loop over nuclides and accumulate fractions for each isotope
    for nuclide in self.nuclides:
      if isinstance(nuclide, Isotope):
        if self.frac_type == 'atom':
          fracs[nuclide.Z] += self.nuclide_atom_fractions[nuclide.zaid]
        elif self.frac_type == 'weight':
          fracs[nuclide.Z] += self.nuclide_weight_fractions[nuclide.zaid]
        rm_zaids.append(nuclide.zaid)
      else:
        new_nuclide_list.append(nuclide)

    # Loop over each isotope and create an element and append it to the new nuclide list
    for Z, frac in fracs.items():
      ZAID = get_elemental_zaid(Z)
      if self.frac_type == 'atom':
        new_nuclide_list.append(Element(Z, atom_fractions={zaid_to_ZA(k)[1] : v for k, v in self.nuclide_atom_fractions.items() if zaid_to_ZA(k)[0] == Z}))
        self.nuclide_atom_fractions[ZAID] = frac
      elif self.frac_type == 'weight':
        new_nuclide_list.append(Element(Z, weight_fractions={zaid_to_ZA(k)[1] : v for k, v in self.nuclide_weight_fractions.items() if zaid_to_ZA(k)[0] == Z}))
        self.nuclide_weight_fractions[ZAID] = frac

    # Delete isotopes from fraction dictionaries
    for zaid in rm_zaids:
      del self.nuclide_atom_fractions[zaid]
      del self.nuclide_weight_fractions[zaid]

    # Set new nuclide list and normalize fractions
    self.nuclides = new_nuclide_list
    self.normalize()

  def normalize(self):
    # Sort the list of nuclides by ZAID
    self.nuclides = sorted(self.nuclides, key=lambda x : int(x.zaid))

    # Calculate the missing atom or weight fraction info
    if self.frac_type == 'atom':
      for n in self.nuclides:
        if n.zaid not in self.nuclide_weight_fractions:
          self.nuclide_weight_fractions[n.zaid] = af_to_wf(n.zaid, self.nuclide_atom_fractions, {k.zaid : k.get_mass() for k in self.nuclides})
    elif self.frac_type == 'weight':
      for n in self.nuclides:
        if n.zaid not in self.nuclide_atom_fractions:
          self.nuclide_atom_fractions[n.zaid] = wf_to_af(n.zaid, self.nuclide_weight_fractions, {k.zaid : k.get_mass() for k in self.nuclides})

    # Normalize the sum of fractions to unity
    atom_frac_sum = sum([self.nuclide_atom_fractions[n.zaid] for n in self.nuclides])
    for n in self.nuclides:
      self.nuclide_atom_fractions[n.zaid] /= atom_frac_sum

    weight_frac_sum = sum([self.nuclide_weight_fractions[n.zaid] for n in self.nuclides])
    for n in self.nuclides:
      self.nuclide_weight_fractions[n.zaid] /= weight_frac_sum

  def to_MCNP_material_card(self, frac_type):
    s = ''
    if self.name is not None:
      s += f'C {self.name}'
    else:
      s += f'C Material {self.number}'
    if self.density is not None:
      s += f' | Density: '
      if self.density > 0:
        s += f'{self.density} [a/(b-cm)]{linesep}'
      else:
        s += f'{abs(self.density)} [g/cc]{linesep}'
    else:
      s += linesep
    if frac_type.lower() == 'atom':
      s += f'{"M" + f"{self.number}":<5s}' + f'{linesep}     '.join((f'{n.zaid:>5s}  {self.nuclide_atom_fractions[n.zaid]:.6E} $ {n.name}' for n in self.nuclides))
    elif frac_type.lower() == 'weight':
      s += f'{"M" + f"{self.number}":<5s}' + f'{linesep}     '.join((f'{n.zaid:>5s}  {-1 * self.nuclide_weight_fractions[n.zaid]:.6E} $ {n.name}' for n in self.nuclides))
    return s

  def get_molar_mass(self) -> float:
    return sum([self.nuclide_atom_fractions[nuclide.zaid] * nuclide.get_mass() for nuclide in self.nuclides])

  def print_table(self):
    mat_dict = {
      'ZAID'          : [n.zaid for n in self.nuclides],
      'Isotope'       : [n.name for n in self.nuclides],
      'Atom Fraction' : [f'{self.nuclide_atom_fractions[n.zaid]:.6E}' for n in self.nuclides],
      'Mass Fraction' : [f'{self.nuclide_weight_fractions[n.zaid]:.6E}' for n in self.nuclides]
    }

    atom_densities = None
    if self.density is not None:
      if self.density > 0:
        atom_densities = {n.zaid : self.density * self.nuclide_atom_fractions[n.zaid] for n in self.nuclides}
        mat_dict['Atom Density [a/(b-cm)]'] = [f'{atom_densities[n.zaid]:.6E}' for n in self.nuclides]
      else:
        atom_densities = {n.zaid : (6.02214076E-1 * abs(self.density) / self.get_molar_mass()) * self.nuclide_atom_fractions[n.zaid] for n in self.nuclides}
        mat_dict['Atom Density [a/(b-cm)]'] = [f'{atom_densities[n.zaid]:.6E}' for n in self.nuclides]

    print(PrintedTable(mat_dict, title=self.name if self.name is not None else '', alignment=['<'] * (5 if self.density is not None else 4)))

    if (self.density is not None) and (atom_densities is not None):
      print(f'Total atom density: {sum(atom_densities.values()):6E} [a/(b-cm)]')

  def write_latex_table(self, file_path, column_specifiers):
    table = {}
    for char in column_specifiers:
      if char == 'Z':
        table['ZAID'] = [n.zaid for n in self.nuclides]
      elif char == 'N':
        table['Isotope' if all(['-' in n.name for n in self.nuclides]) else 'Element'] = [n.name for n in self.nuclides]
      elif char[0] == 'A':
        fmt = char.split(':')[1] if ':' in char else '.6E'
        table['Atom Fraction'] = [to_latex_scientific_notation(f'{self.nuclide_atom_fractions[n.zaid]:{fmt}}') if 'e' in fmt.lower() else f'{self.nuclide_atom_fractions[n.zaid]:{fmt}}' for n in self.nuclides]
      elif char[0] == 'W':
        fmt = char.split(':')[1] if ':' in char else '.6E'
        table['Weight Fraction'] = [to_latex_scientific_notation(f'{self.nuclide_weight_fractions[n.zaid]:{fmt}}') if 'e' in fmt.lower() else f'{self.nuclide_atom_fractions[n.zaid]:{fmt}}' for n in self.nuclides]
    PrintedTable(table).write_latex_longtable(file_path, bold_headers=True)

def mix_materials(mats, fractions, frac_type, num=1):
  new_mat = Material(frac_type=frac_type, number=num)
  for mat, frac in zip(mats, fractions):
    for nuclide in mat.nuclides:
      if frac_type == 'atom':
        new_mat.add_nuclide(nuclide, frac * mat.nuclide_atom_fractions[nuclide.zaid])
      elif frac_type == 'weight':
        new_mat.add_nuclide(nuclide, frac * mat.nuclide_weight_fractions[nuclide.zaid])
  new_mat.normalize()
  return new_mat

def water(num):
  m = Material('atom', number=num, nuclides={'1000' : 2/3, '8000' : 1/3})
  m.split_elements()
  return m

def uo2(num):
  m = Material('atom', number=num, nuclides={'8000' : 2/3, '92234' : 0.00009, '92235' : 0.010124, '92236' : 0.000046, '92238' : 0.323072})
  m.split_elements()
  return m

def graphite_rg(num):
  m = Material('atom', number=num, nuclides={'5000' : 0.000001, '6000' : 0.999999})
  m.split_elements()
  return m

material_fns = {'water' : water, 'uo2' : uo2, 'graphite_rg' : graphite_rg}

def get_compendium_material(library, index: int=None, name: str=None) -> Material:
  """
  Obtain material reference from PNNL's Compendium of Material Composition Data for Radiation Transport Modeling (rev. 2) [PNNL-15870]
  """
  with open(pjoin(dirname(__file__), '..', 'data', f'{library}.yaml'), 'r') as f:
    mat_ref = safe_load(f)
  if (name is not None) and (index is not None):
    raise ValueError('Name and index cannot both be provided!')
  elif (name is None) and (index is None):
    raise ValueError('Must provide either index or name keyword arguments!')
  this_mat = None
  if index is not None:
    this_mat = mat_ref[index]
  elif name is not None:
    for k in mat_ref:
      if mat_ref[k]['name'] == name:
        this_mat = mat_ref[k]
        break
  if this_mat is None:
    raise ValueError('Could not find referenced material!')
  return Material('atom', nuclides=this_mat['composition'], density=-this_mat['mass density'], name=this_mat['name'])
