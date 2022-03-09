from cmath import sqrt


def get_denominator(C, delta_sqrt):
    if abs(C + delta_sqrt) > abs(C - delta_sqrt):
        return C + delta_sqrt
    else:
        return C - delta_sqrt


def parasolve(f, init, eps=1e-3, N=100):
    n = 3
    xs = init
    print("xs: ", xs)
    ys = [f(x) for x in xs]
    print("ys: ", ys)
    while n < N and abs(ys[2]) > eps:
        n += 1
        div_diff_2_1 = (ys[2] - ys[1])/(xs[2] - xs[1])
        div_diff_2_0 = (ys[2] - ys[0])/(xs[2] - xs[0])
        div_diff_1_0 = (ys[1] - ys[0])/(xs[1] - xs[0])
        div_diff_2_1_0 = (div_diff_2_1 - div_diff_1_0)/(xs[2] - xs[0])
        C = div_diff_2_1 + div_diff_2_0 - div_diff_1_0
        delta_sqrt = sqrt(C*C - 4*ys[2]*div_diff_2_1_0)
        denominator = get_denominator(C, delta_sqrt)
        x_n = xs[2] - (2*ys[2]/denominator)
        y_n = f(x_n)
        xs = (xs[1], xs[2], x_n)
        ys = (ys[1], ys[2], y_n)

    if abs(ys[2]) > eps:
        print("maximum number of iterations reached")

    return xs[2]
