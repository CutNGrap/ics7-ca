from prettytable import PrettyTable
from part_1 import *
from part_2 import *


N_MAX = 15
M_MAX = 15


# по y интегрируем гауссом, по x симпсоном
def part_1():
    # x^2+y^2=x <=> (x-0.5)^2+y^2=0.5^2
    # G_integ(sqrt(x^2+y^2)dxdy) = -1_integ_1( 0_integ_2( sqrt(x^2+y^2)dx )dy )
    x_left = 0
    x_right = 2
    y_bottom = -1
    y_top = 1
    table = PrettyTable(['n\m'] + list(range(2, M_MAX)))
    print(f'{gauss_integ(n, m, x_left, x_right, y_bottom, y_top, f):5.5f}')
        table.add_row(string)
    print(table)
    # n = 5
    # m = 10
    # print(gauss_integ(n, m, x_left, x_right, y_bottom, y_top, f))


def part_2():
    table = PrettyTable(['x', 'y', 1, 2, 3, 4, 5])
    table_values = []
    x = [1, 2, 3, 4, 5, 6]
    h = x[1] - x[0]
    y = [0.571, 0.889, 1.091, 1.231, 1.333, 1.412]
    for i in range(len(x)):
        table_values.append([x[i], y[i], 0, 0, 0, 0, 0])
    column_1 = col_1(h, y)
    column_2 = col_2(h, y)
    column_3 = col_3(x, y, column_1)
    column_4 = col_4(x, y)
    column_5 = col_5(h, y)
    columns = [column_1, column_2, column_3, column_4, column_5]
    for j in range(len(columns)):
        for k in range(len(columns[j])):
            table_values[k][2 + j] = f'{columns[j][k]:.5}'

    for i in range(len(x)):
        table.add_row(table_values[i])
    print(table)


part_1()
part_2()
