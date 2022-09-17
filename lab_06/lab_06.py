from part_1 import *
from part_2 import *
from prettytable import PrettyTable

NMAX = 20
MMAX = 20

x = [1, 2, 3, 4, 5, 6]
h = x[1] - x[0]
y = [0.571, 0.889, 1.091, 1.231, 1.333, 1.412]

def integ_1():
    # x^2+y^2=x <=> (x-0.5)^2+y^2=0.5^2
    table = PrettyTable(['n', 'm', 'value'])
    x_left = 0
    x_right = 2
    y_bottom = -1
    y_top = 1
    n = 11
    m =  12
    table.add_row([n, m, f'{double_integr(n, m, x_left, x_right, y_bottom, y_top, f):.5f}'])
    print(table)

def diff_2():
    table = PrettyTable(['x', 'y', 1, 2, 3, 4, 5])
    table_values = []
    for i in range(len(x)):
        table_values.append([x[i], y[i], 0, 0, 0, 0, 0])
    print('\t1 - правая разностная производная')
    print('\t2 - средняя разностная производная')
    print('\t3 - 2ая формула Рунге с правосторонней производной')
    print('\t4 - выравнивающие коэфициенты')
    print('\t5 - вторая разностная производная')
    column_1 = col_1(h, y)
    column_2 = col_2(h, y)
    column_3 = col_3(x, y, column_1)
    column_4 = col_4(x, y)
    column_5 = col_5(h, y)
    columns = [column_1, column_2, column_3, column_4, column_5]
    for j in range(len(columns)):
        for k in range(len(columns[j])):
            table_values[k][2 + j] = '{:.5f}'.format(columns[j][k])

    for i in range(len(x)):
        table.add_row(table_values[i])
    print(table)

if __name__ == "__main__":
    integ_1()
    diff_2()