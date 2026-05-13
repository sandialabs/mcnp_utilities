from mcnp_utilities.lib.tables import PrintedTable, to_latex_scientific_notation

data = {'a' : [1.3, 2.4e-8, 3.7e5], 'b' : ['D', 'E', 'F']}
fd = {'a' : to_latex_scientific_notation, 'b' : lambda x : x.lower()}
PrintedTable(data, column_formatter_dict=fd).write_latex_table('table.tex')
