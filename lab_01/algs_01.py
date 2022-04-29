import math
import copy

EPS = 1e-6

def print_newt_table(table):
    n = len(table)
    for i in range(n):
        for j in range(n - i + 1):
            print("{:<8.3f}".format(table[i][j]), end = ' ')
        print()

def load_dots_from_file(filename):
    try:
        f = open(filename)
    except:
        print("File doesn't exist")
        return []
    dots = []
    line = f.readline()
    while line:
        try:
            x, y, dy = map(float, line.split())
            dots.append([x, y, dy])
        except:
            print("File wasn't properly read")
            break
        line = f.readline()
    f.close()
    return dots

def compare_and_print(y_n, tab_n, y_h, tab_h):
    print("\nComparison table:\n")
    print("{:^15s}{:^15s}".format("Newton", "Hermite"))
    n = len(tab_n)
    m = len(tab_h)
    for i in range(len(tab_h)):
        if i < n - 1:
            cur_n = "{:^8.5f}".format(tab_n[0][i])
        else:
            cur_n = "-"
        print("{:^15s}{:^15.5f}".format(cur_n, tab_h[0][i]))
    print("\ny_newton  = ", y_n, "\ny_hermite = ", y_h)
    print()
    
def is_monotone(dot_arr):
    monotone = 1
    if dot_arr[0][1] <= dot_arr[1][1]:
            is_rising = 1
    else:
        is_rising = -1
    for i in range(len(dot_arr) - 1):
        if (is_rising * (dot_arr[i][1] - dot_arr[i + 1][1]) > 0):
            monotone = 0
            break
    return monotone

def dichotomy(dot_arr, n, key):
    global EPS
    a, b = dot_arr[0][0], dot_arr[len(dot_arr) - 1][0]
    y_a, y_b = dot_arr[0][1], dot_arr[len(dot_arr) - 1][1]
    center = (b + a) / 2
    y_c, t = newton_interpol(dot_arr, center, n) if key != "hermite" else hermite_interpole(dot_arr, center, n)
    while (abs(y_c) > EPS):
        if (y_c * y_a >= 0):
            a = center
        else:
            b = center
        center = (b + a) / 2
        y_c, t = newton_interpol(dot_arr, center, n) if key != "hermite" else hermite_interpole(dot_arr, center, n)
    return center


# ====================================================== == Ньютоновская интерполяция ================================================= ===================
def choose_dots(x, dot_list, n):
    if n > len(dot_list):
        return []
    else:
        before = n // 2
        after = n - before
        arr = []
        last = 0
        for i in range(len(dot_list)):
            if dot_list[i][0] < x:
                last = i
            else:
                break


        if (last + 1) < before:
            before = last + 1
            after = n - before           
        elif (len(dot_list) - last - 1) < after:
            after = len(dot_list) - last - 1
            before = n - after

        for i in dot_list[last - before + 1 : last + after + 1]:
            arr.append(i)
        return arr

def make_newton_table(dot_arr):
    n = len(dot_arr)
    m = len(dot_arr) + 1
    table  = [0] * n
    for i in range(n):
        table[i] = [0] * m
    for i in range(n):
        table[i][0] = dot_arr[i][0]
        table[i][1] = dot_arr[i][1]
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

def newton_interpol(input_list, x, n):
    sorted_list = sorted(input_list, key = lambda x: x[0])
    arr = choose_dots(x, sorted_list, n)
    if len(arr) != n:
        return None, None
    newton_table = make_newton_table(arr)

    return newton_polinome(newton_table, x), newton_table

def reverse_newton(dot_list, n):
    sorted_list = sorted(dot_list, key = lambda x: x[0])
    if is_monotone(sorted_list):
        for i in range(len(sorted_list)):
            t = sorted_list[i][0]
            sorted_list[i][0] = sorted_list[i][1]
            sorted_list[i][1] = t
        sorted_list = sorted(sorted_list, key = lambda x:x[0])
        arr = choose_dots(0, sorted_list, n)
        if len(arr) != n:
            return None
        table = make_newton_table(arr)

        return newton_polinome(table, 0)
    
    else:
        return dichotomy(dot_list, n, "newton")
    

# ====================================================== == Эрмитовская интерполяция ================================================= ===================

def hermite_interpole(input_list, x, n):
    sorted_list = sorted(input_list, key = lambda x: x[0])
    arr = choose_dots(x, sorted_list, n)
    if len(arr) != n:
        return None, None
    hermite_table = make_hermite_table(arr)

    return newton_polinome(hermite_table, x), hermite_table

def make_hermite_table(dot_arr):
    n = len(dot_arr) * 2
    m = n + 2
    table  = [0] * n
    for i in range(n):
        table[i] = [0] * m
    for i in range(n - 1):
        table[i][0] = dot_arr[i//2][0]
        table[i][1] = dot_arr[i//2][1]
        if i % 2 == 0:
            table[i][2] = dot_arr[i//2][2]
        else:
            x1, x2 = dot_arr[i//2][0], dot_arr[i//2 + 1][0]
            y1, y2 = dot_arr[i//2][1], dot_arr[i//2 + 1][1]
            table[i][2] = (y2 - y1) / (x2 - x1)
    table[n - 1][0] = dot_arr[len(dot_arr) - 1][0]
    table[n - 1][1] = dot_arr[len(dot_arr) - 1][1]
    for j in range(3, n + 1):
        for i in range(n - j + 1):
            table[i][j] = (table[i][j - 1] - table[i + 1][j - 1]) / (table[i][0] - table[i + (j - 1)][0])
    return table

def reverse_hermite(dot_list, n):
    sorted_list = sorted(dot_list, key = lambda x: x[0])
    if is_monotone(sorted_list):
        for i in range(len(sorted_list)):
            t = sorted_list[i][0]
            sorted_list[i][0] = sorted_list[i][1]
            sorted_list[i][1] = t
            sorted_list[i][2] = 1e6 if sorted_list[i][2] == 0 else 1 / sorted_list[i][2]
        sorted_list = sorted(sorted_list, key = lambda x:x[0])
        arr = choose_dots(0, sorted_list, n)
        if len(arr) != n:
            return None
        table = make_hermite_table(arr)

        return newton_polinome(table, 0)
    
    else:
        return dichotomy(dot_list, n, "hermite")
 

def main():
    print("Введите название файла, значение x для интерполяции и количество точек для интерполяции через пробел:")
    try:
        filename, x, count = input().split()
        x = float(x)
        count = int(count)
    except:
        print("Неверный ввод")
        return
    input_list = load_dots_from_file(filename)

    y1, newton_table = newton_interpol(input_list, x, count)
    y2, hermite_table = hermite_interpole(input_list, x, count)

    print_table(newton_table, "newton")
    print_table(hermite_table, "hermite")

    compare_and_print(y1, newton_table, y2, hermite_table)

    reverse_list = copy.deepcopy(input_list)
    print("Reverse Newton interpolation:  x = ", reverse_newton(reverse_list, count))
    reverse_list = copy.deepcopy(input_list)
    print("Reverse Hermite interpolation: x = ", reverse_hermite(input_list, count))
    

if __name__ == '__main__':
    main()