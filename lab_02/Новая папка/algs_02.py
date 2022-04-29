import math
import copy
import numpy as np

EPS = 1e-6

def print_table(table, nz, nx, ny):
    for i in range(nz + 1):
        for j in range(ny + 1):
            for k in range(nx + 1):
                print("{:8.4}".format(table[i][j][k]), end = ' ')
            print()
        print("\n")

def print_newt_table(table):
    n = len(table)
    for i in range(n):
        for j in range(n - i + 1):
            print("{:<8.3f}".format(table[i][j]), end = ' ')
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
    return table

def get_last_ind(arr, value):
    last = 0
    for i in range(len(arr)):
        if arr[i] < value:
            last = i
        else:
            break
    return last

def make_n_table(input_tab, y_vals, z_vals, nx, ny, nz, x_mean, y_mean, z_mean):
    table = np.zeros([nz + 1, ny + 1, nx + 1])
    linesy = np.zeros([nz + 1, ny + 1])
    linesz = np.zeros(nz + 1)
    for i in range(nz + 1):
        for j in range(ny + 1):
            index = get_last_ind(input_tab[i][j], x_mean)
            arr = choose_dots(input_tab[i][j], index, nx + 1)

            for k in range(nx + 1):
                table[i][j][k] = arr[k]
    for i in range(nz + 1):
        index = get_last_ind(y_vals[i], y_mean)
        arr = choose_dots(y_vals[i], index, ny + 1)
        for j in range(ny + 1):
            linesy[i][j] = arr[j]

    index = get_last_ind(z_vals, z_mean)
    arr = choose_dots(z_vals, index, nz + 1)
    for i in range(nz + 1):
        linesz[i] = arr[i]

    return table, linesy, linesz

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
def choose_dots(dot_list, index, n):
    arr = []
    before = n // 2
    after = len(dot_list) - before
    if (index + 1) < before:
        before = index + 1
        after = n - before           
    elif (len(dot_list) - index - 1) < after:
        after = len(dot_list) - index - 1
        before = n - after

    for i in dot_list[index - before + 1 : index + after + 1]:
        arr.append(i)
    return arr

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
    for i in range(1, len(table[0]) - 1):
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
        return
    x_vals, y_vals, z_vals = [],[],[]
    table = load_data(filename, x_vals, y_vals, z_vals)
    # print(x_vals, y_vals, z_vals, sep = '\n\n')

    table, y_vals, z_vals = make_n_table(table, y_vals, z_vals, nx, ny, nz, x_mean, y_mean, z_mean)
    x_vals = y_vals[0]
    # print(x_vals, y_vals, z_vals, sep = '\n\n')
    # print_table(table, nz, nx, ny)
    ans = interpol(table, x_vals, y_vals, z_vals, nx, ny, nz, x_mean, y_mean, z_mean)
    # print(x_vals, y_vals, z_vals, sep = '\n\n')
    print(ans)
    

if __name__ == '__main__':
    main()