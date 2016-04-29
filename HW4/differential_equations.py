import math
import matplotlib.pyplot as plt

def euler_generator(f, t0, x0, h):
    """Approximate a first order ODE using Euler's Method.

    This generator produces a sequence of approximations to the
    differential equation dx/dt = f(t, x(t)).  t is the variable, x is
    a function of t and is the solution to the differential equation.

    Args:
        f: the function of t and x(t).
        t0: the initial value of the variable.
        x0: the value of x(t0).
        h: the step size.

    Returns:
        A generator for the approximations.  Each approximation in the
        sequence is separated by the step size, h.
    """
    x = x0
    t = t0
    yield x
    while True:
        x = x + (h*f(t, x))
        t = t + h
        yield x

def mod_euler_generator(f, t0, x0, h):
    """Approximate a first order ODE using the Modified Euler's Method.

    This generator produces a sequence of approximations to the
    differential equation dx/dt = f(t, x(t)).  t is the variable, x is
    a function of t and is the solution to the differential equation.

    Args:
        f: the function of t and x(t).
        t0: the initial value of the variable.
        x0: the value of x(t0).
        h: the step size.

    Returns:
        A generator for the approximations.  Each approximation in the
        sequence is separated by the step size, h.
    """
    x = x0
    t = t0
    yield x
    while True:
        m1 = f(t, x)
        t = t + h
        m2 = f(t, x + (h*m1))
        x = x + (h*0.5*(m1 + m2))
        yield x

def rk4_generator(f, t0, x0, h):
    """Approximate a first order ODE using Runge-Kutta 4.

    This generator produces a sequence of approximations to the
    differential equation dx/dt = f(t, x(t)).  t is the variable, x is
    a function of t and is the solution to the differential equation.

    Args:
        f: the function of t and x(t).
        t0: the initial value of the variable.
        x0: the value of x(t0).
        h: the step size.

    Returns:
        A generator for the approximations.  Each approximation in the
        sequence is separated by the step size, h.
    """
    t = t0
    x = x0
    yield x
    while True:
        k1 = h*f(t, x)
        k2 = h*f(t + (0.5*h), x + (0.5*k1))
        k3 = h*f(t + (0.5*h), x + (0.5*k2))
        k4 = h*f(t + h, x + k3)
        t = t + h
        x = x + (1.0/6)*(k1 + (2*k2) + (2*k3) + k4)
        yield x

def euler_n2_generator(f1, f2, t0, x0, y0, h):
    """Approximate a system of two first order ODE with Euler's Method.

    This generator produces a sequence of approximations to the
    system of differential equations:

        dx/dt = f1(t, x(t), y(t))
        dy/dt = f2(t, x(t), y(t))

    with the initial values t0, x0 = x(t0), y0 = y(t0).

    Args:
        f1: the function equivalent to dx/dt.
        f2: the function equivalent to dy/dt.
        t0: the initial value of the variable.
        x0: the value of x(t0).
        y0: the value of y(t0).
        h: the step size.

    Returns:
        A generator for the approximations.  Each approximation in the
        sequence is separated by the step size, h.  The elements of the
        sequence are tuples (xk, yk), where xk is the approximation for
        x(tk), and similarly for yk.
    """
    t = t0
    x = x0
    y = y0
    yield (x, y)
    while True:
        x = x + (h*f1(t, x, y))
        y = y + (h*f2(t, x, y))
        t = t + h
        yield (x, y)

def rk4_n2_generator(f1, f2, t0, x0, y0, h):
    """Approximate a system of two first order ODE using Runge-Kutta 4.

    This generator produces a sequence of approximations to the
    differential equation dx/dt = f(t, x(t)).  t is the variable, x is
    a function of t and is the solution to the differential equation.

    Args:
        f: the function of t and x(t).
        t0: the initial value of the variable.
        x0: the value of x(t0).
        h: the step size.

    Returns:
        A generator for the approximations.  Each approximation in the
        sequence is separated by the step size, h.
    """
    t = t0
    x = x0
    y = y0
    yield (x, y)
    while True:
        k1 = h*f1(t, x, y)
        m1 = h*f2(t, x, y)
        k2 = h*f1(t + (0.5*h), x + (0.5*k1), y + (0.5*m1))
        m2 = h*f2(t + (0.5*h), x + (0.5*k1), y + (0.5*m1))
        k3 = h*f1(t + (0.5*h), x + (0.5*k2), y + (0.5*m2))
        m3 = h*f2(t + (0.5*h), x + (0.5*k2), y + (0.5*m2))
        k4 = h*f1(t + h, x + k3, y + m3)
        m4 = h*f2(t + h, x + k3, y + m3)
        x = x + (1.0/6)*(k1 + (2.0*k2) + (2.0*k3) + k4)
        y = y + (1.0/6)*(m1 + (2.0*m2) + (2.0*m3) + m4)
        t = t + h
        yield (x, y)

def newtons_method(p0, f, df, tol, max_iter):
    p_old = p0;
    p_new = p0;
    for k in range(max_iter):
        p_new = p_old - (f(p_old) / df(p_old))
        if abs(p_new - p_old) < tol:
            return [True, p_new]
        p_old = p_new
    return [False, p_new]


# Test both Eulers and RK4 using examples 6.4, 6.5, and 6.6 from the
# lecture notes.
example6_4 = iter(euler_generator(lambda t,x: t - x, 0.0, 0.5, 0.2))
example6_5 = iter(mod_euler_generator(lambda t,x: t - x, 0.0, 0.5, 0.2))
example6_6 = iter(rk4_generator(lambda t,x: t - x, 0.0, 0.5, 0.2))
for k in range(6):
    print next(example6_4)
print
for k in range(6):
    print next(example6_5)
print
for k in range(6):
    print next(example6_6)

# Test second order Euler and RK4 using examples from Burden and Faires.
h = 0.1
# From table 5.19 in Burden and Faires ed. 9, chapter 5.9.
#f1 = lambda t,x,y: -4*x + 3*y + 6
#f2 = lambda t,x,y: -2.4*x + 1.6*y + 3.6
#p2 = iter(euler_n2_generator(f1, f2, 0.0, 0.0, 0.0, h))
#p2_1 = iter(rk4_n2_generator(f1, f2, 0.0, 0.0, 0.0, h))

# From table 5.20 in Burden and Faires ed. 9, chapter 5.9.
f1 = lambda t,x,y: y
f2 = lambda t,x,y: (math.e**(2 * t) * math.sin(t)) - (2*x) + (2*y)
p2 = iter(euler_n2_generator(f1, f2, 0.0, -0.4, -0.6, h))
p2_1 = iter(rk4_n2_generator(f1, f2, 0.0, -0.4, -0.6, h))
#f1 = lambda t,x,y: y
#f2 = lambda t,x,y: -3*y - 2*x
#p2 = iter(euler_n2_generator(f1, f2, 0.0, 1.0, 0.0, h))
#p2_1 = iter(rk4_n2_generator(f1, f2, 0.0, 1.0, 0.0, h))
print
for k in range(6):
    print next(p2), next(p2_1)

f1 = lambda t,x,y: y
f2 = lambda t,x,y: -4*(math.pi**2)*x
p2 = iter(euler_n2_generator(f1, f2, 0.0, 0.0, -1.0, h))
p2_1 = iter(rk4_n2_generator(f1, f2, 0.0, 0.0, -1.0, h))
print
for k in range(6):
    print '****'
    print next(p2)
    print next(p2_1)


def plot(x_vals, y_vals, x2_vals, y2_vals, x_min, x_max, y_min, y_max, title):
    plt.axvline(x=0,color='black')
    plt.axhline(y=0,color='black')
    plt.axis([x_min, x_max, y_min, y_max], 1/6)
    plt.plot(x_vals, y_vals, color='red', label='Actual values')
    plt.plot(x2_vals, y2_vals, color='blue', label='Actual values')
    plt.title(title)
    plt.show()

def problem_1():
    f = lambda t,x: t - (2*x)
    t_vals_exact = [k/100.0 for k in range(201)]
    x_vals_exact = [1.25*math.e**(-2*t) + 0.25*(2*t - 1) for t in t_vals_exact]
    #f = lambda t,x: t - x
    #x_vals_exact = [t - 1 + 1.5*(math.e**(-t)) for t in t_vals_exact]
    t0 = 0.0
    x0 = 1.0
    for h in (0.2, 0.1, 0.05):
        title1 = 'Euler with step size ' + str(h)
        title2 = 'Modified Euler with step size ' + str(h)
        title3 = 'RK4 with step size ' + str(h)
        titles = (title1, title2, title3)
        euler = iter(euler_generator(f, t0, x0, h))
        mod_euler = iter(mod_euler_generator(f, t0, x0, h))
        rk4 = iter(rk4_generator(f, t0, x0, h))
        t_vals = [k*h for k in range(int(math.ceil(3.0/h)))]
        x_vals_exact_2 = [1.25*math.e**(-2*t) + 0.25*(2*t - 1) for t in t_vals]
        x_vals = []
        abs_err = []
        iterators = [euler, mod_euler, rk4]
        for sequence, title in zip(iterators, titles):
            for x_exact in x_vals_exact_2:
                x_k = next(sequence)
                x_vals.append(x_k)
                abs_err.append(abs(x_exact - x_k))
            plot(t_vals, x_vals, t_vals_exact, x_vals_exact,
                 0.0, 2.0, 0.0, 1.5, title)
            plot(t_vals, abs_err, [0], [0],
                 0.0, 2.0, 0.0, .12, title + ' absolute error')
            print title, str(max(abs_err))
            x_vals[:] = []
            abs_err[:] = []
        print

def problem_2():
    f1 = lambda t,x,y: y
    f2 = lambda t,x,y: -4 * (math.pi**2) * x
    t_vals_exact = [k/100.0 for k in range(121)]
    x_vals_exact = [-math.cos(2*math.pi*t) for t in t_vals_exact]
    for h in (0.1, 0.05, 0.01, -1):
        title = ''
        if h < 0:
            h = 0.1
            #p2 = iter(rk4_n2_generator(f1, f2, 0.0, 0.0, -1.0, h))
            p2 = iter(rk4_n2_generator(f1, f2, 0.0, -1.0, 0.0, h))
            title = 'RK4 with step size 0.1 absolute error'
        else:
            #p2 = iter(euler_n2_generator(f1, f2, 0.0, 0.0, -1.0, h))
            p2 = iter(euler_n2_generator(f1, f2, 0.0, -1.0, 0.0, h))
            title = 'Euler with step size ' + str(h) + ' absolute error'
        t_vals = [k*h for k in range(int(math.ceil(1.8/h)))]
        x_vals_exact_2 = [-math.cos(2*math.pi*t) for t in t_vals]
        x_vals = []
        abs_err = []
        k = 0
        for x_exact in x_vals_exact_2:
            x_k = list(next(p2))[0]
            x_vals.append(x_k)
            abs_err.append(abs(x_exact - x_k))
            if ((k*h*10) % 1) == 0.0:
                print k*h, x_exact, x_k, abs_err[-1], title
            k = k + 1
        #plot(t_vals, x_vals, t_vals_exact, x_vals_exact,
             #0.0, 1.2, -1.1, 1.1, title)
        #plot(t_vals, abs_err, [0], [0],
             #0.0, 1.2, 0.0, 1e-2, title)

def problem_3():
    f = lambda x: math.sin(x) - (1 + 2*(math.e**(-x)))**(-1)
    df = lambda x: math.cos(x) - (2*(math.e**(-x)))*(1 + 2*(math.e**(-x)))**(-2)
    print
    print newtons_method(0.5, f, df, 1e-11, 10)[1]
    print newtons_method(2.0, f, df, 1e-11, 10)[1]
    f1 = [f(x * 0.01) for x in range(500)]
    x_vals = [x*0.01 for x in range(500)]
    plot(x_vals, f1, [0],  [0], 0.0, 5.0, -2.0, 1.0, 'f(x)')


#problem_1() # TODO: NOTHING!
#problem_2() # TODO: NOTHING!
#problem_3() # TODO: NOTHING!
#problem_4() # TODO: NOTHING!
