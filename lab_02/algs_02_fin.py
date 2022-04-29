import math
import copy
import numpy as np

EPS = 1e-6

def print_table(table, nz, nx, ny):
    for i in range(nz):
        for j in range(ny):
            for k in range(nx):
                print("{:8.4}".format(table[i][j][k]), end = ' ')
            print()
        print("\n")

def print_newt_table(table):
    n = len(table)
    for i in range(n):
        for j in range(n - i + 1):
            print("{:<8.3f}".format(table[i][j]), end = ' ')
        print()
    print()


def load_data(filename, x_vals, y_vals, z_vals):
    try:
        f = open(filename)
    except:
        print("Файл не существует")
        return 
    table_lines = [line.strip() for line in f]
    z = 0
    ammount_y = 5
    cur_ys = []
    nnum = []
    table = []
    for line in table_lines:
        if line == "":
            table.append(nnum)
            nnum = []
        elif "z=" in line:
            z = float(line[list(line).index("=") + 1])
            z_vals.append(z)
        elif "y" in line:
            x_list = list(map(float, line[3:].split()))
            x_vals.append(x_list)
        else:
            linesnum = list(map(float, line.split()))
            nnum.append(linesnum)
            y = linesnum[0]
            cur_ys.append(y)
            if (len(cur_ys) == ammount_y):
                y_vals.append(cur_ys)
                cur_ys = []
            linesnum.pop(0)
    table.append(nnum)
    f.close()
    return table, x_vals[0], y_vals[0], z_vals

def get_last_ind(arr, value):
    last = 0
    for i in range(len(arr)):
        if arr[i] < value:
            last = i
        else:
            break
    return last

def interpol(table, x_vals, y_vals, z_vals, nx, ny, nz, x_mean, y_mean, z_mean):
    yInter = np.zeros([nz + 1, ny + 1])
    for i in range(0, nz + 1):
        for j in range(0 , ny + 1):
            newt_tab = make_newton_table(x_vals, y_vals[i], nx + 1)
            print_newt_table(newt_tab)
            yInter[i][j] = newton_polinome(newt_tab, x_mean)
    print(yInter)
    zInter = []
    for i in range(nz + 1):
        table = make_newton_table(y_vals[i], yInter[i], ny + 1)
        zInter.append(newton_polinome(table, y_mean))
    print(zInter)

    table = make_newton_table(z_vals, zInter, nz + 1)
    result = newton_polinome(table, z_mean)
    return result

# ====================================================== == Ньютоновская интерполяция ================================================= ===================
def choose_dots(dot_x, dot_y, index, n):
    arr_x = []
    arr_y = []
    before = n // 2
    after = n - before
    if (index + 1) < before:
        before = index + 1
        after = n - before           
    elif (len(dot_x) - index - 1) < after:
        after = len(dot_x) - index - 1
        before = n - after

    for i in dot_x[index - before + 1 : index + after + 1]:
        arr_x.append(i)
    for i in dot_y[index - before + 1 : index + after + 1]:
        arr_y.append(i)
    return arr_x, arr_y

def make_newton_table(x, y, n):
    m = len(x) + 1
    table  = [0] * n
    for i in range(n):
        table[i] = [0] * m
    for i in range(n):
        table[i][0] = x[i]
        table[i][1] = y[i]
    for j in range(2, n + 1):
        for i in range(n - j + 1):
            table[i][j] = (table[i][j - 1] - table[i + 1][j - 1]) / (table[i][0] - table[i + (j - 1)][0])
    return table

def newton_polinome(table, x):
    y = 0
    multiplier = 1
    for i in range(1, len(table[0])):
        y += multiplier * table[0][i] 
        multiplier *= x - table[i - 1][0]
    return y

def main():
    print("Введите название файла, значения 3 степени полиномов и значение по x, y, z:")
    try:
        filename, nx, ny, nz, x_mean, y_mean, z_mean = input().split()
        x_mean, y_mean, z_mean = map(float, [x_mean, y_mean, z_mean])
        nx, ny, nz = map(int, [nx, ny, nz])
    except:
        print("Неверный ввод")
        return -1
    x_vals, y_vals, z_vals = [],[],[]
    table, x_vals, y_vals, z_vals = load_data(filename, x_vals, y_vals, z_vals)
    if table == None:
        return -2
    ind_x = get_last_ind(x_vals, x_mean)
    ind_y = get_last_ind(y_vals, y_mean)
    ind_z = get_last_ind(z_vals, z_mean)

    y_interpol = []
    for i in range(len(z_vals)):
        x_interpol = []
        for j in range(len(y_vals)):
            arr_x, arr_y = choose_dots(x_vals, table[i][j], ind_x, nx + 1)
            newt_tab = make_newton_table(arr_x, arr_y, nx + 1)
            x_interpol.append(newton_polinome(newt_tab, x_mean))
        arr_y, arr_x_int = choose_dots(y_vals, x_interpol, ind_y, ny + 1)
        newt_tab = make_newton_table(arr_y, arr_x_int, ny + 1)
        y_interpol.append(newton_polinome(newt_tab, y_mean))
    
    arr_z, arr_y_int = choose_dots(z_vals, y_interpol, ind_z, nz + 1)
    newt_tab = make_newton_table(arr_z, arr_y_int, nz + 1)
    res = newton_polinome(newt_tab, z_mean)
    print(res)
    

if __name__ == '__main__':
    main()