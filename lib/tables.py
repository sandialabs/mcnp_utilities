#!/usr/bin/env python3

from os import linesep, system, remove
from os.path import exists
from pandas import DataFrame as PandasDataFrame
from polars import DataFrame as PolarsDataFrame
from uncertainties import ufloat
from typing import Union, Any


alignment_map = {
  '^' : 'c',
  '<' : 'l',
  '>'  : 'r'
}

special_chars = {
  '\\' : r'\textbackslash',
  '&'  : r'\&',
  '%'  : r'\%',
  '$'  : r'\$',
  '#'  : r'\#',
  '_'  : r'\_',
  '{'  : r'\{',
  '}'  : r'\}',
  '~'  : r'\textasciitilde',
  '^'  : r'\textasciicircum'
}

escape_special_chars_opts = ['none', 'headers', 'elements', 'all']

def to_ufloat_str(num) -> str:
  return fr'${num:.1uS}$'

def to_latex_scientific_notation(num) -> str:
  if isinstance(num, float):
    abcissa, exponent = f'{num:E}'.lower().split('e')
    abcissa = abcissa.strip('0')
    if int(exponent) == 0:
      exponent = '0'
    else:
      exponent = exponent.lstrip('+').lstrip('0') if int(exponent) > 0 else f'-{exponent.lstrip("-").lstrip("0")}'
  elif isinstance(num, str):
    if 'e' in num.lower():
      abcissa, exponent = num.lower().split('e')
      if exponent.startswith('+'):
        exponent = exponent[1:].lstrip('0')
      elif exponent.startswith('-'):
        exponent = f'-{exponent[1:].lstrip("0")}'
    else:
      return num
  return fr'${abcissa} \cdot 10^{{{exponent}}}$'

def to_latex_math(num_str) -> str:
  return fr'${num_str}$'

class PrintedTable:
  def __init__(
      self,
      df:             Union[dict, PandasDataFrame, PolarsDataFrame],
      title:          str         = '',
      alignment:      list[str]   = [],
      center_headers: bool        = False,
      column_formatter_dict: dict = {}
    ) -> None:
    self.title = title
    self.center_headers = center_headers
    if isinstance(df, dict):
      self.df = PandasDataFrame(df)
    else:
      self.df = df
    self.column_formatter_dict = column_formatter_dict
    self.reset_table()
    self.alignment = ['^' for _ in range(len(self.headers))] if not alignment else alignment

  def compute_column_widths(self):
    self.column_widths = [max([len(self.headers[j])] + [len(self.elements[i][j]) for i in range(self.nrows)]) for j in range(self.ncols)]

  def transpose(self, new_first_header=None):
    old_headers = list(self.df.columns)
    new_headers = [new_first_header if new_first_header is not None else self.headers[0]] + list(self.df[self.headers[0]])
    self.df = self.df.transpose()
    self.df.insert(0, new_first_header, old_headers)
    self.df.drop(old_headers[0], inplace=True)
    self.df.columns = [str(h) for h in new_headers]
    self.column_formatter_dict = {}
    self.reset_table()
    self.alignment = ['^' for _ in range(len(self.headers))]

  def reset_table(self):
    self.headers = [str(h) for h in self.df.columns]
    self.column_formatter_dict = {c : str if c not in self.column_formatter_dict else self.column_formatter_dict[c] for c in self.headers}
    self.nrows, self.ncols = len(self.df), len(self.headers)
    self.elements = [[self.column_formatter_dict[c](self.df.iloc[i][c] if isinstance(self.df, PandasDataFrame) else self.df[c][i]) for c in self.headers] for i in range(len(self.df))]
    self.compute_column_widths()

  def bold_headers(self):
    self.headers = [f'\\textbf{{{h}}}' for h in self.headers]
    self.compute_column_widths()

  def escape_special_chars(self, which='headers'):
    for special_char, replacement in special_chars.items():
      if (which == 'all') or (which == 'headers'):
        for j in range(self.ncols):
          self.headers[j] = self.headers[j].replace(special_char, replacement) if special_char in self.headers[j] else self.headers[j]
      if (which == 'all') or (which == 'elements'):
        for i in range(self.nrows):
          self.elements[i] = [ele.replace(special_char, replacement) if special_char in ele else ele for ele in self.elements[i]]
    self.compute_column_widths()

  def elements_in_math_mode(self):
    self.elements = [[f'${ele}$' for ele in self.elements[i]] for i in range(self.nrows)]
    self.compute_column_widths()

  def set_column(self, column: list[Any], index: int):
    if len(self.elements) != len(column):
      raise AssertionError(f'Column length ({len(column)}) must be the same as the number of rows in the table ({len(self.elements)})!')
    for i in range(len(self.elements)):
      self.elements[i][index] = column[i]
    self.compute_column_widths()

  def get_column(self, index):
    return [self.elements[i][index] for i in range(len(self.elements))]

  def to_latex_table(
      self,
      label:                   str  = '',
      indentation_level:       int  = 0,
      indentation_size:        int  = 2,
      tabular_env:             str  = 'tabular',
      tabular_env_opts:        str  = '',
      label_indentifier:       str  = 'tab',
      table_placement_options: str  = '!htbp',
      table_type:              str  = 'builtin',
      hlines:                  bool = True,
      outer_box:               bool = False,
      bold_headers:            bool = False,
      stretch:                 str  = '',
      resize:                  str  = None,
      standalone:              bool = False,
      escape_special_chars:    str  = 'headers',
      elements_in_math_mode:   bool = False
    ) -> str:

    if table_type not in ['builtin', 'booktabs']:
      raise AssertionError(f'Unrecognized table type "{table_type}"!')

    if escape_special_chars not in escape_special_chars_opts:
      raise AssertionError(f'Argument "escape_special_chars" must be one of {escape_special_chars_opts}!')

    if escape_special_chars != 'none':
      self.escape_special_chars(escape_special_chars)
    if elements_in_math_mode:
      self.elements_in_math_mode()
    if bold_headers:
      self.bold_headers()
    bar_str, lines = r" \hline" if hlines else "", []

    indent_level = 0

    if standalone:
      lines.append(r'\documentclass{standalone}')
      lines.append(r'\begin{document}')
    else:
      lines.append(f'\\begin{{table}}{f"[{table_placement_options}]" if table_placement_options else ""}')
      indent_level += 1
      if stretch:
        lines.append(rf'{" " * indentation_size * indent_level}\renewcommand{{\arraystretch}}{{{stretch}}}')
      lines.append(rf'{" " * indentation_size * indent_level}\centering')
      if self.title:
        lines.append(fr'{" " * indentation_size * indent_level}\caption{{{self.title}}}')
      if label:
        lines.append(rf'{" " * indentation_size * indent_level}\label{{{label_indentifier}:{label}}}')
    if resize is not None:
      lines.append(f'{" " * indentation_size * indent_level}\\resizebox{{{resize}\\textwidth}}{{!}}{{%')
      indent_level += 1
    lines.append(f'{" " * indentation_size * indent_level}\\begin{{{tabular_env}}}{"" if not tabular_env_opts else f"[{tabular_env_opts}]"}{{{("|" if outer_box else "") + "|".join([alignment_map[al] for al in self.alignment]) + ("|" if outer_box else "")}}}')
    indent_level += 1
    if table_type == 'booktabs':
      lines.append(f'{" " * indentation_size * indent_level}\\toprule')
    elif outer_box:
      lines.append(f'{" " * indentation_size * indent_level}\\hline')
    lines.append(f'{" " * indentation_size * indent_level}{" & ".join(f"{self.headers[j]:<{self.column_widths[j]}}" for j in range(self.ncols))} \\\\{bar_str if not table_type == "booktabs" else ""}')
    if table_type == 'booktabs':
      lines.append(f'{" " * indentation_size * indent_level}\\midrule')
    for i in range(self.nrows):
      lines.append(f'{" " * indentation_size * indent_level}{" & ".join(f"{self.elements[i][j]:<{self.column_widths[j]}}" for j in range(self.ncols))} \\\\{bar_str if i+1 != self.nrows else ""}')
    if table_type == 'booktabs':
      lines.append(f'{" " * indentation_size * indent_level}\\bottomrule')
    elif outer_box:
      lines.append(f'{" " * indentation_size * indent_level}\\hline')
    indent_level -= 1
    lines.append(f'{" " * indentation_size * indent_level}\\end{{{tabular_env}}}{"%" if resize is not None else ""}')
    indent_level -= 1
    if resize is not None:
      lines.append(f'{" " * indentation_size * indent_level}}}')
    if standalone:
      lines.append(r'\end{document}')
    else:
      lines.append(rf'\end{{table}}')

    self.reset_table()
    return linesep.join([f'{(" " * indentation_size) * indentation_level}{line}' for line in lines])

  def to_latex_subtable(
      self,
      width:                 str  = r'\textwidth',
      label:                 str  = '',
      tabular_env:           str  = 'tabular',
      tabular_env_opts:      str  = '',
      indentation_level:     int  = 0,
      indentation_size:      int  = 2,
      label_indentifier:     str  = 'stab',
      table_type:            str  = 'builtin',
      hlines:                bool = True,
      outer_box:             bool = False,
      bold_headers:          bool = False,
      stretch:               str  = '',
      escape_special_chars:  str  = 'headers',
      elements_in_math_mode: bool = False
    ) -> str:

    if table_type not in ['builtin', 'booktabs']:
      raise AssertionError(f'Unrecognized table type "{table_type}"!')

    if escape_special_chars not in escape_special_chars_opts:
      raise AssertionError(f'Argument "escape_special_chars" must be one of {escape_special_chars_opts}!')

    if escape_special_chars != 'none':
      self.escape_special_chars(escape_special_chars)
    if elements_in_math_mode:
      self.elements_in_math_mode()
    if bold_headers:
      self.bold_headers()
    bar_str, lines = r" \hline" if hlines else "", []

    lines.append(f'\\begin{{subtable}}[*]{{{width}}}')
    if stretch:
      lines.append(rf'  \renewcommand{{\arraystretch}}{{{stretch}}}')
    lines.append(rf'  \centering')
    if self.title:
      lines.append(fr'  \caption{{{self.title}}}')
    if label:
      lines.append(rf'  \label{{{label_indentifier}:{label}}}')
    lines.append(f'  \\begin{{{tabular_env}}}{"" if not tabular_env_opts else f"[{tabular_env_opts}]"}{{{("|" if outer_box else "") + "|".join([alignment_map[al] for al in self.alignment]) + ("|" if outer_box else "")}}}')
    if table_type == 'booktabs':
      lines.append(r'    \toprule')
    elif outer_box:
      lines.append(r'    \hline')
    lines.append(f'    {" & ".join(f"{self.headers[j]:<{self.column_widths[j]}}" for j in range(self.ncols))} \\\\{bar_str if not table_type == "booktabs" else ""}')
    if table_type == 'booktabs':
      lines.append(r'    \midrule')
    for i in range(self.nrows):
      lines.append(f'    {" & ".join(f"{self.elements[i][j]:<{self.column_widths[j]}}" for j in range(self.ncols))} \\\\{bar_str if i+1 != self.nrows else ""}')
    if table_type == 'booktabs':
      lines.append(r'    \bottomrule')
    elif outer_box:
      lines.append(r'    \hline')
    lines.append(f'  \\end{{{tabular_env}}}')
    lines.append(rf'\end{{subtable}}')

    self.reset_table()
    return linesep.join([f'{(" " * indentation_size) * indentation_level}{line}' for line in lines])

  def to_latex_longtable(
      self,
      label:                   str  = '',
      indentation_level:       int  = 0,
      indentation_size:        int  = 2,
      label_indentifier:       str  = 'tab',
      table_placement_options: str  = 'c',
      table_type:              str  = 'builtin',
      hlines:                  bool = True,
      outer_box:               bool = False,
      bold_headers:            bool = False,
      copy_headers:            bool = False,
      stretch:                 str  = '',
      escape_special_chars:    str  = 'headers',
      elements_in_math_mode:   bool = False
    ) -> str:

    if table_type not in ['builtin', 'booktabs']:
      raise AssertionError(f'Unrecognized table type "{table_type}"!')

    if escape_special_chars not in escape_special_chars_opts:
      raise AssertionError(f'Argument "escape_special_chars" must be one of {escape_special_chars_opts}!')

    if escape_special_chars != 'none':
      self.escape_special_chars(escape_special_chars)
    if elements_in_math_mode:
      self.elements_in_math_mode()
    if bold_headers:
      self.bold_headers()
    bar_str, lines = r" \hline" if hlines else "", []

    lines.append(f'\\begin{{longtable}}{f"[{table_placement_options}]" if table_placement_options else ""}{{{("|" if outer_box else "") + "|".join([alignment_map[al] for al in self.alignment]) + ("|" if outer_box else "")}}}')
    if stretch:
      lines.append(rf'  \renewcommand{{\arraystretch}}{{{stretch}}}')
    lines.append(rf'  \centering')
    if self.title:
      lines.append(fr'  \caption{{{self.title}}}')
    if label:
      lines.append(rf'  \label{{{label_indentifier}:{label}}}')
    if table_type == 'booktabs':
      lines.append(r'  \toprule')
    elif outer_box:
      lines.append(r'  \hline')
    lines.append(f'  {" & ".join(f"{self.headers[j]:<{self.column_widths[j]}}" for j in range(self.ncols))} \\\\{bar_str if not table_type == "booktabs" else ""}')
    if table_type == 'booktabs':
      lines.append(r'  \midrule')
    if copy_headers:
      lines.append(r'  \endfirsthead')
      if table_type == 'booktabs':
        lines.append(r'  \toprule')
      elif outer_box:
        lines.append(r'  \hline')
      lines.append(f'  {" & ".join(f"{self.headers[j]:<{self.column_widths[j]}}" for j in range(self.ncols))} \\\\{bar_str if not table_type == "booktabs" else ""}')
      if table_type == 'booktabs':
        lines.append(r'  \midrule')
      lines.append(r'  \endhead')
    for i in range(self.nrows):
      lines.append(f'  {" & ".join(f"{self.elements[i][j]:<{self.column_widths[j]}}" for j in range(self.ncols))} \\\\{bar_str if i+1 != self.nrows else ""}')
    if table_type == 'booktabs':
      lines.append(r'  \bottomrule')
    elif outer_box:
      lines.append(r'  \hline')
    lines.append(rf'\end{{longtable}}')

    self.reset_table()
    return linesep.join([f'{(" " * indentation_size) * indentation_level}{line}' for line in lines])

  def write_latex_table(self, out_file, **kwargs) -> None:
    with open(out_file, 'w') as f:
      f.write(self.to_latex_table(**kwargs))

  def write_latex_subtable(self, out_file, **kwargs) -> None:
    with open(out_file, 'w') as f:
      f.write(self.to_latex_subtable(**kwargs))

  def write_latex_longtable(self, out_file, **kwargs) -> None:
    with open(out_file, 'w') as f:
      f.write(self.to_latex_longtable(**kwargs))

  def create_pdf_table(self, out_file_base, **kwargs) -> None:
    self.write_latex_table(f'{out_file_base}.tex', standalone=True, **kwargs)
    system(f'pdflatex {out_file_base}.tex > /dev/null')
    for fyle in [f'{out_file_base}.fls', f'{out_file_base}.fdb_latexmk', f'{out_file_base}.log', f'{out_file_base}.aux', f'{out_file_base}.tex']:
      if exists(fyle):
        remove(fyle)

  def __str__(self) -> str:
    lines = []
    # Create table header
    table_header = '│'
    for i, header in enumerate(self.headers):
      table_header += f' {header:{self.alignment[i] if not self.center_headers else "^"}{self.column_widths[i]}} ┆'
    table_header = table_header[:-1] + '│'
    table_width = len(table_header)

    divider = '┌'
    for v in self.column_widths:
      divider += ('─' * (v + 2)) + '┬'
    divider = divider[:-1] + '┐'
      # Print title
    if self.title:
      lines.append(f'{self.title:^{table_width}}')
    lines.append(divider)
    # Print the header
    lines.append(table_header)
    divider = '╞'
    for v in self.column_widths:
      divider += ('═' * (v + 2)) + '╪'
    divider = divider[:-1] + '╡'
    lines.append(divider)
    # Print the body
    for i in range(self.nrows):
      line = '│'
      for j in range(self.ncols):
        line += f' {self.elements[i][j]:{self.alignment[j]}{self.column_widths[j]}} ┆'
      line = line[:-1] + '│'
      lines.append(line)
    divider = '└'
    for v in self.column_widths:
      divider += ('─' * (v + 2)) + '┴'
    divider = divider[:-1] + '┘'
    lines.append(divider)

    return linesep.join(lines)

def write_latex_table(
    tables:                  list[PrintedTable],
    out_file:                str,
    stab_kwargs:             Union[list[dict], dict] = [],
    indentation_level:       int                     = 0,
    indentation_size:        int                     = 2,
    title:                   str                     = '',
    label:                   str                     = '',
    table_placement_options: str                     = '!htbp',
    stretch:                 str                     = ''
  ) -> None:
  lines = []
  lines.append(f'\\begin{{table}}{f"[{table_placement_options}]" if table_placement_options else ""}')
  if stretch:
    lines.append(rf'  \renewcommand{{\arraystretch}}{{{stretch}}}')
  lines.append(rf'  \centering')
  if title:
    lines.append(fr'  \caption{{{title}}}')
  if label:
    lines.append(rf'  \label{{tab:{label}}}')
  for i, table in enumerate(tables):
    if stab_kwargs:
      if isinstance(stab_kwargs, list):
        lines.append(table.to_latex_subtable(**stab_kwargs[i], indentation_level=indentation_level+1, indentation_size=indentation_size))
      elif isinstance(stab_kwargs, dict):
        lines.append(table.to_latex_subtable(**stab_kwargs, indentation_level=indentation_level+1, indentation_size=indentation_size))
    else:
      lines.append(table.to_latex_subtable(indentation_level=indentation_level+1, indentation_size=indentation_size))
  lines.append(rf'\end{{table}}')
  with open(out_file, 'w') as f:
    [f.write(line) for line in linesep.join([f'{(" " * indentation_size) * indentation_level}{line}' for line in lines])]
