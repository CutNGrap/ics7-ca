from math import *
import matplotlib.pyplot
import pandas as pd
import numpy as np
import random

SIZE = 50
DELTA = 1e-7
FILE = "points.csv"
FILE_2DIM = "points_2dim.csv"

def one_dimensional():
    print('''
    1 - без учёта весов
    2 - с учётом весов
    
    0 - отмена
    ''')
    decision = int(input('Ваш выбор: '))
    if decision == 0:
        return
    n = int(input('Введите степень аппроксимирующего полинома: '))
    x, y, p = loadTableXY()
    N = len(x)
    gauss_matrix = generate_gauss_matrix_one_dim(x, y, p, n, N, decision == 1)
    a = solve_gauss_matrix(gauss_matrix)
    approx_x = [x[0]]
    dx = x[-1] - x[0]
    step = dx / (len(x) * 3)
    for i in range(len(x) * 3):
        approx_x.append(approx_x[-1] + step)
    approx_y = []
    for i in range(len(approx_x)):
        approx_y.append(0)
        for j in range(len(a)):
            approx_y[-1] += a[j] * approx_x[i] ** j
    x_med = medium(x, p, decision == 1)
    y_med = medium(y, p, decision == 1)
    matplotlib.pyplot.plot([x_med], [y_med], "r*")
    matplotlib.pyplot.plot(x, y, "*", color="green")
    matplotlib.pyplot.plot(approx_x, approx_y)
    matplotlib.pyplot.show()

def read_two_dim_point():
    data = pd.read_csv(FILE_2DIM)
    data.head(data.size)
    x = data['x'].values
    y = data['y'].values
    z = data['z'].values
    r = data["r"].values
    return x, y, z, r


def two_dimensional():
    print('''
        1 - без учёта весов
        2 - с учётом весов

        0 - отмена
        ''')
    decision = int(input('Ваш выбор: '))
    if decision == 0:
        return
    n = int(input('Введите степень аппроксимирующего полинома(поддерживается только 1 или 2): '))
    x, y, z, p = read_two_dim_point()
    N = len(x)
    X, Y = get_grid(x, y)
    if n == 1:
        gauss_matrix = get_gauss_two_dim_first(x, y, z, p, N, decision == 1)
        coef = solve_gauss_matrix(gauss_matrix)
        Z = coef[0] + coef[1] * X + coef[2] * Y
    elif n == 2:
        gauss_matrix = get_gauss_two_dim_second(x, y, z, p, N, decision == 1)
        coef = solve_gauss_matrix(gauss_matrix)
        Z = coef[0] + coef[1] * X + coef[2] * Y + coef[3] * X * Y + coef[4] * X * X + coef[5] * Y * Y
        print(coef)
    else:
        return
    ax = matplotlib.pyplot.axes(projection='3d')
    ax.plot_surface(X, Y, Z, color='red')
    ax.scatter(x, y, z, color='blue')
    matplotlib.pyplot.show()

def get_grid(x, y):
    approx_x = [x[0]]
    dx = x[-1] - x[0]
    step = dx / (len(x) * 3)
    for i in range(len(x) * 3):
        approx_x.append(approx_x[-1] + step)
    sorted(y)
    approx_y = [y[0]]
    dy = y[-1] - y[0]
    step = dy / (len(y) * 3)
    for i in range(len(y) * 3):
        approx_y.append(approx_y[-1] + step)
    return np.meshgrid(x, y)

def generate_gauss_matrix_one_dim(x, y, p, n, N, same_weight):
    gauss_matrix = []
    for m in range(n + 1):
        gauss_string = []
        for k in range(n + 1):
            summa = 0
            for i in range(N):
                pi = 1 if same_weight else p[i]
                summa += pi * x[i] ** (k + m)
            gauss_string.append(summa)
        summa = 0
        for i in range(N):
            pi = 1 if same_weight else p[i]
            summa += pi * y[i] * x[i] ** m
        gauss_string.append(summa)
        gauss_matrix.append(gauss_string)
    return gauss_matrix

def solve_gauss_matrix(matrix):
    for i in range(len(matrix) - 1):
        if matrix[i][i] == 0:
            k = i + 1
            while k < len(matrix) and matrix[k][i] != 0:
                k += 1
            if k == len(matrix):
                return
            else:
                swap_lines(matrix, i, k)
        for j in range(i + 1, len(matrix)):
            # if abs(matrix[j][i]) < DELTA:
                coef = -(matrix[i][i] / matrix[j][i])
                mult_line(matrix[j], coef)
                matrix[j] = add_lines(matrix[j], matrix[i])
    x = [0] * len(matrix)
    for i in range(len(matrix) - 1, -1, -1):
        summa = matrix[i][-1]
        for j in range(i + 1, len(matrix)):
            summa -= matrix[i][j] * x[j]
        x[i] = summa / matrix[i][i]
    return x


def get_gauss_two_dim_first(x, y, z, p, N, same_weight):
    gauss_matrix = []
    for i in range(3):
        gauss_matrix.append([0, 0, 0, 0])
    for i in range(N):
        pi = 1 if same_weight else p[i]
        values = [1, x[i], y[i], z[i]]
        diffs = [1, x[i], y[i]]
        for j in range(3):
            for k in range(4):
                gauss_matrix[j][k] += pi * values[k] * diffs[j]
    return gauss_matrix


def get_gauss_two_dim_second(x, y, z, p, N, same_weight):
    gauss_matrix = []
    for i in range(6):
        gauss_matrix.append([0, 0, 0, 0, 0, 0, 0])
    for i in range(N):
        pi = 1 if same_weight else p[i]
        values = [1, x[i], y[i], x[i] * y[i], x[i] * x[i], y[i] * y[i], z[i]]
        diffs = [1, x[i], y[i], x[i] * y[i], x[i] * x[i], y[i] * y[i]]
        for j in range(6):
            for k in range(7):
                gauss_matrix[j][k] += pi * values[k] * diffs[j]
    return gauss_matrix


def mult_line(line, coef):
    for i in range(len(line)):
        line[i] *= coef

def add_lines(line1, line2):
    result = [0] * len(line1)
    for i in range(len(line1)):
        result[i] = line1[i] + line2[i]
    return result

def medium(x, p, same_weight):
    chis = 0
    znam = 0
    for i in range(len(x)):
        pi = 1 if same_weight else p[i]
        chis += pi * x[i]
        znam += pi
    return chis / znam

def loadTableXY():
    data = pd.read_csv(FILE)
    data.head(data.size)
    x = data['x'].values
    y = data['y'].values
    ro = data["r"].values
    return x, y, ro
    
def create_table():
    global FILE
    print('''
    1 - для одной переменной
    2 - для двух переменных
    
    0 - отмена
    ''')
    fileindex = int(input('Ваш выбор: '))
    if fileindex == 0:
        return
    N = int(input('Введите количество узлов: '))
    if fileindex == 1:
        f = open("randompoint.csv", "w")
        ymax = 6
        ymin = -5
        y0 = 0
        xarr = np.zeros(N)
        yarr = np.zeros(N)
        r = np.zeros(N)
        for i in range(N):
            y = y0 + random.random() * (ymax - ymin) + ymin
            y0 = y
            xarr[i] = i
            yarr[i] = y
            r[i] = random.choice(range(1, 51))
        d = {'x':xarr, 'y':yarr, 'r':r}
        df = pd.DataFrame(data = d)
        df.to_csv(path_or_buf = f)
        return xarr, yarr, r
    elif fileindex == 2:
        f = open("randompoint_2dim.csv", 'w')
        xarr = np.zeros(N)
        yarr = np.zeros(N)
        zarr = np.zeros(N)
        r = np.zeros(N)
        size = SIZE
        for i in range(N):
            xarr[i] = i
            yarr[i] = random.choice(range(size + 1))
            zarr[i] = random.choice(range(size + 1))
            r[i] = random.choice(range(1, 51))
        d = {'x':xarr, 'y':yarr, 'z':zarr, 'r':r}
        df = pd.DataFrame(data = d)
        df.to_csv(path_or_buf = f)
        return xarr, yarr, r

if __name__ == "__main__":
    code = 1
    while code != 0:
        print('''
        1 - сгенерировать таблицу
        2 - провести одономерную аппроксимацию
        3 - провести двумерную аппроксимацию
        
        0 - выход
        ''')
        code = int(input('Ваш выбор: '))
        if code == 1:
            create_table()
        elif code == 2:
            one_dimensional()
        elif code == 3:
            two_dimensional()

