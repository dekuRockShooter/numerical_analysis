import math

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
        x = x + (1.0/6)*(k1 + (2*k2) + (2*k3) + k4)
        y = y + (1.0/6)*(m1 + (2*m2) + (2*m3) + m4)
        t = t + h
        yield (x, y)

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
