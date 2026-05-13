#!/usr/bin/env python3

import unittest
from pandas import DataFrame as PandasDataFrame
from polars import DataFrame as PolarsDataFrame
from mcnp_utilities.lib.tables import PrintedTable


class TestPrintedTable(unittest.TestCase):
  def test_different_init_table(self):
    my_pt = PrintedTable(df=PandasDataFrame({'h1' : ['a', 'b', 'c'], 'h2' : [1, 2, 3]}))
    df_pandas_init = my_pt.to_latex_table()

    my_pt = PrintedTable(df=PolarsDataFrame({'h1' : ['a', 'b', 'c'], 'h2' : [1, 2, 3]}))
    df_polars_init = my_pt.to_latex_table()

    self.assertMultiLineEqual(df_pandas_init, df_polars_init)

  def test_different_init_subtable(self):
    my_pt = PrintedTable(df=PandasDataFrame({'h1' : ['a', 'b', 'c'], 'h2' : [1, 2, 3]}))
    df_pandas_init = my_pt.to_latex_subtable(r'\linewidth')

    my_pt = PrintedTable(df=PolarsDataFrame({'h1' : ['a', 'b', 'c'], 'h2' : [1, 2, 3]}))
    df_polars_init = my_pt.to_latex_subtable(r'\linewidth')

    self.assertMultiLineEqual(df_pandas_init, df_polars_init)

  def test_different_init_longtable(self):
    my_pt = PrintedTable(df=PandasDataFrame({'h1' : ['a', 'b', 'c'], 'h2' : [1, 2, 3]}))
    df_pandas_init = my_pt.to_latex_longtable()

    my_pt = PrintedTable(df=PolarsDataFrame({'h1' : ['a', 'b', 'c'], 'h2' : [1, 2, 3]}))
    df_polars_init = my_pt.to_latex_longtable()

    self.assertMultiLineEqual(df_pandas_init, df_polars_init)

if __name__ == '__main__':
  unittest.main()
