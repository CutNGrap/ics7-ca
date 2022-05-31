import numpy

DELTA = 1e-7


class Polynome:
    def __init__(self, coef):
        self.coef = coef

    def __str__(self):
        string = ''
        for i in range(len(self.coef)):
            coef = f'{self.coef[i]:5.2f}'
            if not self.coef[i] < 0:
                coef = f'+{self.coef[i]:5.2f}'
            if i == len(self.coef) - 1:
                string += f'{coef} '
            else:
                string += f'{coef}*x^{len(self.coef) - i - 1} '
        return string

    def get_coeffs(self):
        return self.coef

    def dp_dy(self):
        if len(self.coef) == 1:
            new_coeffs = [0]
        else:
            new_coeffs = [0] * (len(self.coef) - 1)
        for i in range(len(self.coef) - 1):
            new_coeffs[i] = self.coef[i] * (len(self.coef) - i - 1)
        return Polynome(new_coeffs)

    def roots(self):
        if len(self.coef) == 1:
            return False
        elif len(self.coef) == 2:
            return [self.coef[1] / self.coef[0]] if abs(self.coef[0]) > DELTA else False
        else:
            deriv_pol = self.dp_dy()
            roots = deriv_pol.roots()
            if type(roots) == list:
                roots = [-1] + roots + [1]
                new_roots = [0] * (len(roots) - 1)
                for i in range(len(roots) - 1):
                    xl = roots[i]
                    xr = roots[i + 1]
                    med_x = (xl + xr) / 2
                    med_y = self.value(med_x)
                    fl = True
                    while abs(med_y) > DELTA and fl:
                        if med_y * self.value(xl) > 0:
                            xl = med_x
                        else:
                            xr = med_x
                        med_x = (xl + xr) / 2
                        fl = med_y
                        med_y = self.value(med_x)
                        fl = abs(1 - med_y / fl) > DELTA
                    new_roots[i] = med_x
                return new_roots
            else:
                return False

    def value(self, x):
        y = 0
        for i in range(len(self.coef)):
            y += self.coef[i] * x ** (len(self.coef) - i - 1)
        return y


class Legander_polynome(Polynome):
    def __init__(self, n):
        coefs = self.__generate__(n)
        super().__init__(coefs)

    def __generate__(self, n):
        if n == 0:
            return [1]
        elif n == 1:
            return [1, 0]
        else:
            coefs_n_1 = self.__generate__(n - 1)
            coefs_n_2 = self.__generate__(n - 2)
            for i in range(len(coefs_n_1)):
                coefs_n_1[i] *= 2 - 1 / n
            coefs_n_1.append(0)
            for i in range(len(coefs_n_2)):
                coefs_n_1[2 + i] -= coefs_n_2[i] * (1 - 1 / n)
            return coefs_n_1


def generate_matrix(n, roots):
    a = [[0] * n for i in range(n)]
    b = [0] * n
    for i in range(n):
        for j in range(len(roots)):
            a[i][j] = roots[j] ** i
        if i % 2 == 0:
            b[i] = 2 / (i + 1)
        else:
            b[i] = 0
    return a, b


def gauss_integ(n, m, x_left, x_right, y_bottom, y_top, function):
    pol = Legander_polynome(n)
    roots = pol.roots()
    print(roots)
    a, b = generate_matrix(n, roots)
    a_i = list(numpy.linalg.solve(a, b))
    s = 0
    for i in range(n):
        fun = lambda x: function(x, (y_bottom + y_top) / 2 + (y_top - y_bottom) / 2 * roots[i])
        s += a_i[i] * simpson_integ(m, x_left, x_right, fun)
    s *= (y_top - y_bottom) / 2
    return s


def simpson_integ(n, bot, top, function):
    h = (top - bot) / n
    s = 0
    cur_step = bot
    for i in range(n // 2):
        s += function(cur_step) + 4 * function(cur_step + h) + function(cur_step + 2 * h)
        cur_step += 2 * h
    s *= h / 3
    return s


# вспомогательная функция; sqrt(x^2+y^2) внутри x^2+y^2=2x
def f(x, y):
    if (x - 1) ** 2 + y ** 2 - 1 > DELTA:
        return 0
    else:
        return (x ** 2 + y ** 2) ** 0.5

if __name__ == "__main__":
    x_left = 0
    x_right = 2
    y_bottom = -1
    y_top = 1
    n = 10
    m =  5
    print(gauss_integ(n, m, x_left, x_right, y_bottom, y_top, f))