from matplotlib import pyplot as plt
import numpy as np

EPS = 1e-4
l = 10
T0 = 300
R = 0.5
F0 = -10
a1 = 0.0134
b1 = 1
c1 = 4.35e-4
m1 = 1
alpha0 = 1.94e-2
sigma = 1.5e3
gamma = 0.2e-2
N = 50

h = l / N

def init():
    ev = dict()

    ev['A'] = [0] * (N + 1)
    ev['B'] = [0] * (N + 1)
    ev['C'] = [0] * (N + 1)
    ev['D'] = [0] * (N + 1)

    ev['A1'] = [0] * (N + 1)
    ev['B1'] = [0] * (N + 1)
    ev['C1'] = [0] * (N + 1)
    ev['D1'] = [0] * (N + 1)

    ev['ksi'] = [0] * (N + 1)
    ev['teta'] = [0] * (N + 1)

    ev['k'] = [0] * (N + 1)
    ev['p'] = [0] * (N + 1)
    ev['alpha'] = [0] * (N + 1)
    ev['f'] = [0] * (N + 1)
    ev['T'] = [float(T0)] * (N + 1) 
    return ev

def func_p(t):
    return 2 / R * func_alpha(t)

def func_k(t):
    return a1 * (b1 + c1 * t**m1)

def func_f(t):
    return 2 * T0 / R * func_alpha(t)

def func_alpha(t):
    return alpha0 * (t / sigma - 1) ** 4 + gamma

def func_d_k(t):
    return a1 * c1 * m1 * t**(m1-1)

def func_d_p(t):
    return 2 / R * func_d_alpha(t)

def func_d_f(t):
    return 2 * T0 / R * func_d_alpha(t)

def func_d_alpha(t):
    return (4 * alpha0 * (t / sigma - 1)**3) / sigma

def K_plus_half(k, i):
    return (k[i] + k[i + 1]) / 2

def K_minus_half(k, i):
    return (k[i - 1] + k[i]) / 2



def init_left(ev):
    T = ev['T']
    ev['k'][1] = func_k(T[1])

    ev['A'][0] = 0
    ev['B'][0] = K_plus_half(ev['k'], 0) + ev['p'][0] * h**2
    ev['C'][0] = K_plus_half(ev['k'], 0)
    ev['D'][0] = ev['f'][0] * h**2 + F0 * h

    A_d = 0
    B_d = lambda x: func_d_k(x) / 2 + func_d_p(x) * h**2
    C_d = lambda x: func_d_k(x) / 2
    D_d = lambda x: func_d_f(x) * h**2

    ev['A1'][0] = 0
    ev['B1'][0] = B_d(T[0]) * T[0] + ev['B'][0] - C_d(T[0]) * T[1] - D_d(T[0])
    ev['C1'][0] = -C_d(T[1]) * T[0] + C_d(T[1]) * T[1]  + ev['C'][0]
    ev['D1'][0] = -ev['B'][0] * T[0] + ev['C'][0] * T[1] + ev['D'][0]

def init_right(ev):
    T = ev['T']

    ev['A'][N] = K_minus_half(ev['k'], N)
    ev['B'][N] = K_minus_half(ev['k'], N) + ev['p'][N] * h**2 + ev['alpha'][N] * h
    ev['C'][N] = 0.0
    ev['D'][N] = ev['f'][N] * h**2 + ev['alpha'][N] * T0 * h

    A_d = lambda x: func_d_k(x) / 2
    B_d = lambda x: func_d_k(x) / 2 + func_d_p(x) * h**2 + func_d_alpha(x) * h
    C_d = 0.0
    D_d = lambda x: func_d_f(x) * h**2 + func_d_alpha(x) * T0 * h

    ev['A1'][N] = A_d(T[N - 1])  * T[N - 1] + ev['A'][N] - A_d(T[N - 1]) * T[N]
    ev['B1'][N] = - A_d(T[N]) * T[N - 1] + B_d(T[N]) * T[N] + ev['B'][N] - D_d(T[N])
    ev['C1'][N] = 0.0
    ev['D1'][N] = ev['A'][N] * T[N - 1] - ev['B'][N] * T[N] + ev['D'][N]

def count(ev, x_range):
    max_dy = 1
    T = ev["T"]

    while max_dy > EPS:
        # print(*T)
        init_left(ev)
        for i in range(1, N):
            ev['k'][i] = func_k(T[i])
            ev['alpha'][i] = func_alpha(T[i])
            ev['p'][i] = func_p(T[i])
            ev['f'][i] = func_f(T[i])

            ev['k'][i + 1] = func_k(T[i + 1])

            ev['A'][i] = K_minus_half(ev['k'], i)
            ev['B'][i] = K_minus_half(ev['k'], i) +  K_plus_half(ev['k'], i) + ev['p'][i] * h**2
            ev['C'][i] = K_plus_half(ev['k'], i)
            ev['D'][i] = ev['f'][i] * h**2

            A_d = lambda T: func_d_k(T) / 2
            B_d = lambda T: func_d_k(T) + func_d_p(T) * h ** 2
            C_d = lambda T: func_d_k(T) / 2
            D_d = lambda T: func_d_f(T) * h ** 2

            ev['A1'][i] = A_d(T[i - 1]) * T[i - 1] + ev['A'][i] - A_d(T[i - 1]) * T[i]
            ev['B1'][i] = -A_d(T[i]) * T[i - 1] + B_d(T[i]) * T[i] + ev['B'][i] - C_d(T[i]) * T[i + 1] - D_d(T[i])
            ev['C1'][i] = -C_d(T[i + 1]) * T[i] + C_d(T[i + 1]) * T[i + 1] + ev['C'][i]
            ev['D1'][i] = ev['A'][i] * T[i - 1] - ev['B'][i] * T[i] + ev['C'][i] * T[i + 1] + ev['D'][i]

            ev['ksi'][i] = ev['C1'][i - 1] / (ev['B1'][i - 1] - ev['A1'][i - 1] * ev['ksi'][i - 1])
            ev['teta'][i] = (ev['A1'][i - 1] * ev['teta'][i - 1] + ev['D1'][i - 1]) / (ev['B1'][i - 1] - ev['A1'][i - 1] * ev['ksi'][i - 1])

        init_right(ev)
        i = N
        ev['ksi'][i] = ev['C1'][i - 1] / (ev['B1'][i - 1] - ev['A1'][i - 1] * ev['ksi'][i - 1])
        ev['teta'][i] = (ev['A1'][i - 1] * ev['teta'][i - 1] + ev['D1'][i - 1]) / (ev['B1'][i - 1] - ev['A1'][i - 1] * ev['ksi'][i - 1])


        dy = [0] * (N + 1)
        max_dy = 0

        for i in range(N - 1, -1, -1):
            if i != N:
                dy[i] = ev['ksi'][i + 1] * dy[i + 1] + ev['teta'][i + 1]
            else:
                dy[i] = (ev['A1'][i] * ev['teta'][i] + ev['D1'][i]) / (ev['B1'][i] - ev['A1'][i] * ev['ksi'][i])
            if T[i] != 0:
                max_dy = max(max_dy, abs(dy[i] / T[i]))
            else:
                max_dy = max(max_dy, 0)
        if max_dy > EPS:
            for i in range(N + 1):
                T[i] += dy[i]

if __name__ == "__main__":
    ev = init()
    x_range = np.arange(0, l + h, h)
    count(ev, x_range)

    ax = plt.subplot()

    for i in range(0, N+1, 10):
        print('{:E}'.format(ev['T'][i]))
    ax.plot(x_range, ev['T'], '-')
    ax.set_xlabel('Расстояние х от левого конца стержня, см')
    ax.set_ylabel('Температура, К')
    ax.grid(True)
    plt.show()