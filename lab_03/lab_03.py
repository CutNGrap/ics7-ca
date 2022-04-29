import matplotlib.pyplot as plt

from algs_02_fin import *

FILE = 'data.txt'
VAR = 2


def read_points():
    f = open(FILE, 'r')
    i = 0
    data = []
    line = f.readline()
    while line:
        try:
            point = list(map(float, line.split()))
            data.append([point[0], point[1]])
        except:
            print('Неверные входные данные')
            f.close()
            return []
        i += 1
        line = f.readline()
    f.close()
    return data


def diff_2_newton(data, necessary_point, newt_tab):
    y_0_2 = newt_tab[0][3]
    y_0_3 = newt_tab[0][4]
    return 2 * y_0_2 + (6 * necessary_point - 2 * (data[0][0] + data[1][0]) - 2 * data[2][0]) * y_0_3


def find_h(data, N):
    h = [0] * (N - 1)
    for i in range(1, N):
        h[i - 1] = data[i][0] - data[i - 1][0]
    return h

def find_a(data, N):
    a = [0] * (N - 1)
    for i in range(N - 1):
        a[i] = data[i][1]
    return a

def find_d(N, h, c):
    d = [0] * (N - 1)
    for i in range(N - 1):
        d[i] = (c[i + 1] - c[i]) / (3 * h[i])
    return d


def find_b(data, N, h, c):
    b = [0] * (N - 1)
    for i in range(N - 1):
        b[i] = (data[i + 1][1] - data[i][1]) / h[i] - h[i] * ((c[i + 1] + 2 * c[i]) / 3)
    return b

def find_c(data, N, h, c1 = 0, cN_plus_1 = 0):
    c = [0] * N
    E = [0] * (N - 1)
    n = [0] * (N - 1)

    E[0] = 0
    n[0] = c1
    for i in range(1, N - 1):
        f = 3 * ((data[i + 1][1] - data[i][1]) / h[i] - (data[i][1] - data[i - 1][1]) / h[i - 1])
        n[i] = (f - h[i - 1] * n[i - 1]) / (h[i - 1] * E[i - 1] + 2 * (h[i - 1] + h[i]))
        E[i] = (-h[i]) / (h[i - 1] * E[i - 1] + 2 * (h[i - 1] + h[i]))

    c[-1] = cN_plus_1
    if c[-1]:
        alpha = h[N - 2] * h[N - 3] / 3
        beta = 2 * h[N - 2] * (h[N - 3] + h[N - 2])
        gamma = data[N - 1][1] - data[N - 2][1] + h[N - 2] * (data[N - 3][1] - data[N - 2][1]) / h[N - 3] \
                - c[-1] * h[N - 2] ** 2 / 6
        c[-2] = (gamma - alpha * n[-2]) / (alpha * E[-2] + beta)
    else:
        c[-2] = E[-1] * c[-1] + n[-1]
    for i in range(N - 3, -1, -1):
        c[i] = E[i] * c[i + 1] + n[i]
    return c

def graph(data, N, a, b, c, d):
    x = []
    y = []
    list_x = [0] * 20
    list_y = [0] * 20
    for i in range(N - 1):
        step = (data[i + 1][0] - data[i][0]) / 19
        cur_x = data[i][0]
        for j in range(20):
            list_x[j] = cur_x
            list_y[j] = a[i] + b[i] * (cur_x - data[i][0]) + c[i] * (cur_x - data[i][0]) ** 2 + d[i] * (
                        cur_x - data[i][0]) ** 3
            cur_x += step
        x += list_x
        y += list_y
    return x, y


if __name__ == '__main__':
    data = read_points()
    
    N = len(data)
    h = find_h(data, N)
    a = find_a(data, N)

    if VAR == 0:
        c = find_c(data, N, h)

    elif VAR == 1:
        x = [i[0] for i in data]
        y = [i[1] for i in data]

        arr_x, arr_y = choose_dots(x,y,0, 4)
        newt_tab = make_newton_table(arr_x, arr_y, 4)
        c1 = diff_2_newton(data, data[0][0], newt_tab)
        c = find_c(data, N, h, c1=c1)
    elif VAR == 2:
        x = [i[0] for i in data]
        y = [i[1] for i in data]

        arr_x, arr_y = choose_dots(x,y,0, 4)
        newt_tab = make_newton_table(arr_x, arr_y, 4)
        c1 = diff_2_newton(data, data[0][0], newt_tab)

        arr_x, arr_y = choose_dots(x,y,N, 4)
        newt_tab = make_newton_table(arr_x, arr_y, 4)
        cN_plus_1 = diff_2_newton(data, data[0][0], newt_tab)
        c = find_c(data, N, h, c1=c1, cN_plus_1=cN_plus_1)
    else:
        exit(-1)  

    b = find_b(data, N, h, c)
    d = find_d(N, h, c)
    plt.plot([point[0] for point in data], [point[1] for point in data])
    plt.plot([point[0] for point in data], [point[1] for point in data], "*", color="red")
    x, y = graph(data, N, a, b, c, d)
    plt.plot(x, y)
    plt.show()