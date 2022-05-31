from prettytable import PrettyTable

# правосторонняя разностная производная O(h)
def col_1(h, y):
    dy = [0] * len(y)
    for i in range(len(y) - 1):
        dy[i] = (y[i + 1] - y[i]) / h
    dy[-1] = dy[-2]
    return dy


# центральная разностная производная O(h^2)
def col_2(h, y):
    dy = [0] * len(y)
    dy[0] = (-3 * y[0] + 4 * y[1] - y[2]) / (2 * h)
    for i in range(1, len(y) - 1):
        dy[i] = (y[i + 1] - y[i - 1]) / (2 * h)
    dy[-1] = (3 * y[-1] - 4 * y[-2] + y[-3]) / (2 * h)
    return dy


# 2я формула Рунге с правосторонней производной O(h^2)
def col_3(x, y, column_1):
    m = 2
    p = 1
    h = m * (x[1] - x[0])
    dy_h1 = col_1(h, y[::2])
    dy_h2 = col_1(h, y[1::2])
    dy = [0] * (len(dy_h1) + len(dy_h2))
    for i in range(len(dy_h1) + len(dy_h2)):
        if i % 2 == 0:
            dy[i] = dy_h1[i // 2]
        else:
            dy[i] = dy_h2[i // 2]
    for i in range(len(dy)):
        dy[i] = column_1[i] + (column_1[i] - dy[i]) / (m ** p - 1)
    return dy


# с использованием выравнивающих переменных ksi = y/x ; eta = y; =>
# => y = (a0 * x) / (a1 + a2 * x) => eta = a0/a2 - a1/a2 * ksi
def col_4(x, y):
    eta = [0] * len(y)
    ksi = [0] * len(y)
    for i in range(len(x)):
        eta[i] = y[i]
        ksi[i] = y[i] / x[i]

    deta = [0] * len(eta)
    for i in range(len(eta) - 1):
        deta[i] = (eta[i + 1] - eta[i]) / (ksi[i + 1] - ksi[i])
    deta[-1] = deta[-2]

    dy = [0] * len(y)
    for i in range(len(dy)):
        dy[i] = (y[i] / (x[i] ** 2) * deta[i]) / (deta[i] / x[i] - 1)
    return dy


# вторая разностная производная O(h^2)
def col_5(h, y):
    dy = [0] * len(y)
    dy[0] = (2 * y[0] - 5 * y[1] + 4 * y[2] - y[3]) / (h ** 2)
    for i in range(1, len(y) - 1):
        dy[i] = (y[i - 1] - 2 * y[i] + y[i + 1]) / (h ** 2)
    dy[-1] = (2 * y[-1] - 5 * y[-2] + 4 * y[-3] - y[-4]) / (h ** 2)
    return dy

if __name__ == "__main__":
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
            table_values[k][2 + j] = '{:.5f}'.format(columns[j][k])

    for i in range(len(x)):
        table.add_row(table_values[i])
    print(table)