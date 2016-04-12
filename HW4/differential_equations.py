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

example64 = iter(euler_generator(lambda t,x: t - x, 0.0, 0.5, 0.2))
example65 = iter(mod_euler_generator(lambda t,x: t - x, 0.0, 0.5, 0.2))
example66 = iter(rk4_generator(lambda t,x: t - x, 0.0, 0.5, 0.2))
#example62 = iter(euler_generator(lambda x,t: x - (2*t), 0, 1.0, 0.2))
for k in range(6):
    print next(example64)
print
for k in range(6):
    print next(example65)
print
for k in range(6):
    print next(example66)
