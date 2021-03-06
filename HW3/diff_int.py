import math
import matplotlib.pyplot as plt

def _central_diff(x, f):
    """Return a function that computes the central difference of f.

    The central difference is defined as (f(x + h) - f(x - h)) / 2h.

    Args:
        x: the real number to calculate the derivative of f at.
        f: the function to differentiate.

    Returns:
        g(h): A funtion that calculates the central difference of f(x).
              Its argument, 'h', is the step size used to compute the
              central difference.
    """
    return lambda h: 1.0*(f(x + h) - f(x - h)) / (2*h)
    
def _trapezoid_approx(x0, xn, f):
    """Return a function that approximats a definite integral.

    The function that is returned uses the Trapezoid Rule to evaluate
    a definite integral.

    Args:
        x0: the lower bound of the interval of integration.
        xn: the upper bound of the interval of integration.
        f: the function to integrate.

    Returns:
        g(h): A funtion that approximates the integral of f(x) from x0
              to xn.  The approximation uses the Trapezoid Rule.  The
              function's argument, 'h', is the step size.  It is the
              length of each trapezoid used in the approximation.
              The size of 'h' is assumed to split the interval of
              integration into equal subintervals.
    """
    edge_sum = f(x0) + f(xn)
    def intermediate_sum(h):
        c = 2
        sum = 0
        xk = x0
        while xk < xn - h:
            xk = xk + h
            sum = sum + f(xk)
            c = c + 1
        print c
        return ((sum*2) + edge_sum) * h / 2
    return intermediate_sum

def _approximate(f, h, tol, max_iter):
    """Implement the iterative step in Richardson Extrapolation.

    The step used to get successively better approximations is the
    same for both Richardson Extrapolation and Romberg Integration.
    This function exists in order to remove the code duplication.
    
    It is intended to be used with _central_diff and _trapezoid_approx.
    f should be the function returned by one of those two functions.
    _approximate will then be approximating a derivative at some
    value, or a definite integral between some bounds.

    Args:
        f: the function used in the approximation.  It should take a
           single argument, h, that is the step size.  Calling f(h)
           would be the same as calling f(x + h), where x is given.
           See the example for more detail.
        h: the initial step size.  As the algorithm progresses, this
           is halved, in accordance to the procedure used in Richardson
           extrapolation and Romberg Integration.
        tol: the tolerance for a good approximation.
        max_iter: the maximum number of iterations to perform in case
                  no approximation is good enough.

    Returns:
        The evaluation of f(h) if it is within the tolerance.
        0 is no such value was found.
    """
    M = []
    for k in range(max_iter):
        M.append([0 for _ in range(k + 1)])
        for j in range(k + 1):
            if j == 0:
                M[k][j] = f(h / (2**k))
            else:
                M[k][j] = 1.0*M[k][j-1] + ((M[k-1][j-1] - M[k][j-1]) /
                                           (2**(2*j)-1))
                if abs(M[k][j] - M[k][j - 1]) < tol:
                    return M[k][j]
    return 0

# takes a long time for small h.
def romberg_integration(x0, xn, f, h, tol, max_iter):
    """Evaluate the definate integral of a function.
    
    This uses Romberg Integration to approximate the integral of f(x)
    from x0 to xn.

    Args:
        x0: the lower bound of the interval of integration.
        xn: the upper bound of the interval of integration.
        f: the single variable function to integrate.  It should
           have the form f(x).
        h: the initial step size to use in the approximation.
        tol: the tolerance used to deem an approximation as good.
        max_iter: the maximum number of iterations to perform before
                  aborting the approximation.

    Returns:
        The integral of f(x) from x0 to xn.
        0 if a good approximation was not found.
    """
    trap = _trapezoid_approx(x0, xn, f)
    return _approximate(trap, h, tol, max_iter)
    
def gaussian_quadrature(f, n):
    """Approximates a definite integral using Gaussian quadrature.
    
    The interval of integration is bounded by -1 and 1.  If the integral
    to approximate uses another interval, then a linear transformation to
    convert the interval to [-1, 1] should be performed prior to calling this
    function.
    
    Args:
        f: the integrand defined on the interval [-1, 1].  This should be
           of the form f(x), where x is a real number in [-1, 1].
        n: The number of points to use in the approximation.  The more points,
           the better the approximation.
    
    Returns:
        The approximate integral of f over [-1, 1].
    """
    weights = [
            0.4179591836734693877551020408163265306122448979591836734693877551020408163265306122448979591836734693877551020408163265306122448979591836734693877551020408163265306122448979591836734693877551020408163265306122448979591836734693877551020408163265306122448979591836734693877551020408163265306,
    0.3818300505051189449503697754889751338783650835338627347510834510307055464341297083486846593440448014503146717645853573344928956776383837562443187566373816994263513750309425122069048082192405967657155458166140211350441276674773890998286086501179973611850420964132316867109449892362441359044150514997832,
    0.3818300505051189449503697754889751338783650835338627347510834510307055464341297083486846593440448014503146717645853573344928956776383837562443187566373816994263513750309425122069048082192405967657155458166140211350441276674773890998286086501179973611850420964132316867109449892362441359044150514997832,
    0.2797053914892766679014677714237795824869250652265987645370140326936188104305626768132409429011976187663233752133720515191356369795631199443713526578123368545563592025336909192643194833502493348267909279642862302108156020621692804889667628308214167200738365553989662119973317061306305537859150514997832,
    0.2797053914892766679014677714237795824869250652265987645370140326936188104305626768132409429011976187663233752133720515191356369795631199443713526578123368545563592025336909192643194833502493348267909279642862302108156020621692804889667628308214167200738365553989662119973317061306305537859150514997832,
    0.1294849661688696932706114326790820183285874022599466639772086387246552349720423087156254181629208450894844020016344278810653448938189044626496347079992610378540241163129175889369389737366325173870853629537936262051606784336186365336536081108973206126186723685959653665979209749176350897029150514997832,
    0.1294849661688696932706114326790820183285874022599466639772086387246552349720423087156254181629208450894844020016344278810653448938189044626496347079992610378540241163129175889369389737366325173870853629537936262051606784336186365336536081108973206126186723685959653665979209749176350897029150514997832
    ]
    points = [
                0,
    0.4058451513773971669066064120769614633473820140993701263870432517946638132261256553283126897277465877652867586660480186780142389774087899602458293459431152403705864850136028192946798646997494188869169765542654505357384603100658598476270710450994883480024599267113885472679490162043321422574150514997832,
    -0.4058451513773971669066064120769614633473820140993701263870432517946638132261256553283126897277465877652867586660480186780142389774087899602458293459431152403705864850136028192946798646997494188869169765542654505357384603100658598476270710450994883480024599267113885472679490162043321422574150514997832,
    -0.7415311855993944398638647732807884070741476471413902601199553519674298746721805137928268323668632470596925180931120142436000543982298353471703857152740498332960747607976107150698769026932844561958151246095962171815950287169821619140709720118875391555834601955414971467103462901278094572097150514997832,
    0.7415311855993944398638647732807884070741476471413902601199553519674298746721805137928268323668632470596925180931120142436000543982298353471703857152740498332960747607976107150698769026932844561958151246095962171815950287169821619140709720118875391555834601955414971467103462901278094572093150514997832,
    -0.949107912342758524526189684047851262400770937670617783548769103913063330354840140805730770027925724144300739666995216194195625811353553118277789915859810085013901000179888247732305040104815148851112904940437420579459979108498442397952261081440138823188704950068274774322776063669713039873415051499783203,
    0.949107912342758524526189684047851262400770937670617783548769103913063330354840140805730770027925724144300739666995216194195625811353553118277789915859810085013901000179888247732305040104815148851112904940437420579459979108498442397952261081440138823188704950068274774322776063669713039873415051499783203
    ]
    return sum([w * f(p) for w, p in zip(weights, points)])
    
def plot(x_vals, y_vals, x_min, x_max, y_min, y_max, title):
    plt.axvline(x=0,color='black')
    plt.axhline(y=0,color='black')
    plt.axis([x_min, x_max, y_min, y_max], 1/6)
    plt.plot(x_vals, y_vals, color='red', label='Actual values')
    plt.title(title)
    plt.show()

# Romberg_a,0.318269086201,65
# Romberg_b,0.716371650486,97
# Romberg_c,0.225718812875,129
# Trap_a,0.318309870579,2049
# Trap_b,0.716362807919,1537
# Trap_c,0.225705853263,2049
# Gauss_a,0.318309886184,4
# Gauss_b,0.716362794363,6
# Gauss_c,0.225705832851,7
# Table 1 shows the results of each approximation, accurate to seven decimal
# places.  Using the Trapezoid Rule is far and away the slowest algorithm.
# It, however, agrees with the results from Gaussian quadrature for the firsr
# seven digits after the decimal point.  Romberg integration is comparably much
# faster, but its approximations do not agree as closely to any other algorithm
# the way the Trapezoid Rule and Gaussian Quadrature do.  This is odd since
# Romberg integration should be more accurate than the Trapezoid Rule.  This is
# most likely due to the step sizes used.
def problem_3():
    print '***Problem 3***'
    function_a = lambda x: math.sin(math.pi * x)
    function_a_gauss = lambda u: math.sin(math.pi * (u + 1) * 0.25) * 0.25
    function_b = lambda x: math.pow(math.e, -(x**2.0) / 2.0)
    def function_b_gauss(u):
        return 0.75 * math.pow(math.e, -((0.75*u + 1.25)**2) * 0.5)

    def function_c(x):
        if x == 0:
            return 1
        return math.sin(math.pi * 2.0 * x) / (math.pi * 2.0 * x)

    def function_c_gauss(u):
        if u == -1:
            return 1
        return 0.5 * math.sin(math.pi*(u + 1)) / (math.pi*(u + 1))

    x_vals = [k*0.005 for k in range(101)]
    plot(x_vals, [function_a(x) for x in x_vals], 0.0, 0.5, -0.1, 1.1,
         'sin(pi*x)')
    x_vals = [k*0.005 for k in range(100, 401)]
    plot(x_vals, [function_b(x) for x in x_vals], 0.5, 2.03, -0.1, 1.1,
         'exp(-(x**2) / 2)')
    x_vals = [k*0.005 for k in range(201)]
    plot(x_vals, [function_c(x) for x in x_vals], -0.1, 1.1, -0.3, 1.1,
         'Sinc(2*pi*x)')
    print 'Romberg a: '\
          + str(romberg_integration(0, 0.5, function_a, .01, 1e-7, 30))
    print 'Romberg b: '\
          + str(romberg_integration(0.5, 2, function_b, .01, 1e-7, 30))
    print 'Romberg c: '\
          + str(romberg_integration(0.0, 1.0, function_c, .5, 1e-7, 30))
    function_a_trap = _trapezoid_approx(0.0, 0.5, function_a)
    function_b_trap = _trapezoid_approx(0.5, 2.0, function_b)
    function_c_trap = _trapezoid_approx(0.0, 1.0, function_c)

    def eval(f):
        step = 0.5
        prev_val = f(step)
        diff = 1.0
        while diff > 1e-7:
            print '***'+ str(step)
            step = step / 2.0
            cur_val = f(step)
            diff = abs(cur_val - prev_val)
            prev_val = cur_val
        return cur_val

    print 'Trapezoid a: ' + str(eval(function_a_trap))
    print 'Trapezoid b: ' + str(eval(function_b_trap))
    print 'Trapezoid c: ' + str(eval(function_c_trap))
    print 'Gaussian a: ' + str(gaussian_quadrature(function_a_gauss, 5))
    print 'Gaussian b: ' + str(gaussian_quadrature(function_b_gauss, 5))
    print 'Gaussian c: ' + str(gaussian_quadrature(function_c_gauss, 5))

problem_3()
