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


example62 = iter(euler_generator(lambda x,t: x - (2*t), 0, 1.0, 0.2))
for k in range(6):
    print next(example62)
