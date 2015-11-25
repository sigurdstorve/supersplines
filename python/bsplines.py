import numpy as np




def float_is_zero(v, eps=1e-6):
    """ Test if floating point number is zero. """
    return abs(v) < eps

def special_div(num, den):
    """ Return num/dev with the special rule
    that 0/0 is 0. """
    if float_is_zero(num) and float_is_zero(den):
        return 0.0
    else:
        return num / den

def B(j, p, x, knots):
    """ Compute B-splines using recursive definition. """        
    if p == 0:
        if knots[j] <= x < knots[j+1]:
            return 1.0
        else:
            return 0.0
    else:
        left = special_div((x-knots[j])*B(j,p-1,x,knots), knots[j+p]-knots[j])
        right = special_div((knots[j+1+p]-x)*B(j+1,p-1,x,knots), knots[j+1+p]-knots[j+1])
        return left + right

def get_mu(x, knots, mu=None):
    """
    Find the knot interval index (mu) for a given value x
    This is a prerequisite for alg 2.20 and 2.21.
    Optional to supply a suggestion to try first.
    Throws on error.
    """
    if mu != None:
        if x >= knots[mu] and x < knots[mu+1]:
            return mu
    for mu in range(0, len(knots)-1):
        if x >= knots[mu] and x < knots[mu+1]:
            return mu
    raise RuntimeError("Illegal knot vector or x value")

def R_mat(k, mu, knots, x):
    """
    Compute k'th B-spline matrix (k=1...)
    mu: so that t_{mu} <= x < t_{mu+1}
    knots: the knot vector
    x: value to evaluate in
    """
    res = np.zeros((k,k+1))
    for row in range(k):
        common_denom = knots[mu+(row+1)] - knots[mu+(row+1)-k]
        res[row,row] = special_div(knots[mu+(row+1)] - x, common_denom)
        res[row,row+1] = special_div(x - knots[mu + (row+1) - k], common_denom)
    return res

def alg_220(p, knots, control_points, x, mu=None):
    """
    Compute the spline function value at x
    p: degree
    knots: knots vector
    control_points:
    x: value to evaluate at
    """
    if mu == None:
        # auto-determine mu
        mu = get_mu(x, knots)
    c = np.array(control_points[(mu-p):(mu+1)])
    assert( len(c) == p+1 )
    for k in range(p, 0, -1):   # k = p...1
        c = R_mat(k, mu, knots, x).dot(c)
    return c    

def alg_221(p, knots, x):
    """
    Evaluate all p+1 nonzero B-splines at x
    p: degree
    knots: knot vector
    x: value to evaluate at
    """
    mu = get_mu(x, knots)
    B = np.array([1.0])
    for k in range(1, p+1):
        B = B.dot(R_mat(k, mu, knots, x))
    return B
    
def render_spline(p, knots, control_points, ts):
    """
    Compute points on a spline function using the straightforward
    implementation of the recurrence relation for the B-splines.
    """
    
    ys = []
    for t in ts:
        y = 0.0 
        for j in range(0, len(control_points)):
            y += B(j, p, t, knots)*control_points[j]
        ys.append(y)
    return ys

def render_spline_alg220(p, knots, control_points, ts):
    # TODO: verify input data
    # TODO: reuse get_mu() and alg_220()
    xs = []
    mu = None
    for t in ts:
        # Try the last mu value first, and then search if needed
        if mu == None or not (t >= knots[mu] and t < knots[mu+1]):
            mu = None
            for i in range(0, len(knots)-1):
                if t >= knots[i] and t < knots[i+1]:
                    mu = i
            if mu == None: raise Exception("Unable to determine mu")

        # Compute the B-spline vector
        B_vec = np.array([[1]])
        for k in range(1, p+1):
            Rk = R_mat(k, mu, knots, t)
            B_vec = B_vec.dot(Rk)

        # Dot with correct part of control point vector
        c0 = control_points[(mu-p):(mu+1)]
        x = B_vec.dot(c0)[0]
        xs.append(x)
    return xs

def uniform_regular_knot_vector(n, p, t0=0.0, t1=1.0):
    """
    Create a p+1-regular uniform knot vector for
    a given number of control points
    Throws if n is too small
    """

    # The minimum length of a p+1-regular knot vector
    # is 2*(p+1)
    if n < p+1:
        raise RuntimeError("Too small n for a uniform regular knot vector")

    # p+1 copies of t0 left and p+1 copies of t1 right
    # but one of each in linspace
    return [t0]*p + list(np.linspace(t0, t1, n+1-p)) + [t1]*p

def control_points(p, knots):
    """
    Return the control point abscissa for the control polygon
    of a one-dimensional spline.
    """
    knots = np.array(knots)
    abscissas = []
    for i in range(len(knots)-p-1):
        part = knots[(i+1):(i+1+p)]
        abscissas.append(np.mean(part))
    return abscissas

def render_tensor_prod_spline(p1, p2, knots1, knots2, control_points, us, vs):
    """ Slow ref. impl. """
    res = np.empty((len(us), len(vs)))
    n1, n2 = control_points.shape
    assert(n1 == len(knots1) - p1 - 1)
    assert(n2 == len(knots2) - p2 - 1)
    for u_i,u in enumerate(us):
        for v_i,v in enumerate(vs):
            s = 0.0
            for i in range(n1):
                for j in range(n2):
                    s += control_points[i,j]*B(i, p1, u, knots1)*B(j, p2, v, knots2)
            res[u_i, v_i] = s
    return res

def least_squares_surface_fit(xs, ys, zs, ws1, ws2, p1, p2, knots1, knots2):
    """
    Data points (xs[i], ys[j], zs[i,j]) (i=1..m1, j=1..m2)
    ws1[i] weights (i=1..m1)
    ws2[j] weights (j=1..m2)
    Spline spaces S1 and S2 implicitly defined by 
    (p1, knots1) and (p2, knots2):
        n1=len(knots1)-(p1+1)
        n2=len(knots2)-(p2+1)
    Returns n1 x n2 matrix of coefficients
    """
    m1, m2 = zs.shape
    assert(m1 == len(xs))
    assert(m2 == len(ys))
    assert(len(ws1) == m1)
    assert(len(ws2) == m2)
    
    n1 = len(knots1) - (p1+1)
    n2 = len(knots2) - (p2+1)
    
    # Create matrix A
    A = np.empty((m1, n1))
    for i in range(m1):
        for q in range(n1):            
            A[i, q] = np.sqrt(ws1[i])*B(q, p1, xs[i], knots1)

    # Create matrix _B
    _B = np.empty((m2, n2))
    for j in range(m2):
        for r in range(n2):
            _B[j, r] = np.sqrt(ws2[j])*B(r, p2, ys[j], knots2)
    
    # Create matrix G
    G = np.empty((m1, m2))
    for i in range(m1):
        for j in range(m2):
            G[i, j] = np.sqrt(ws1[i])*np.sqrt(ws2[j])*zs[i, j]

    # Compute the coefficient matrix from eq. (7.20)
    A_trans = A.transpose()
    M1 = np.linalg.inv(A_trans.dot(A))
    M2 = np.linalg.inv(_B.transpose().dot(_B))
    C = M1.dot(A_trans).dot(G).dot(_B).dot(M2)
    
    return C

def lsq_spline_fit(xs, ys, knots, p, ws=None):
    """
    Returns spline coefficients for least-squares
    fit of function samples (x_i, y_i) on an arbitrary
    knot vector and degree.
    """
    if ws == None:
        ws = np.ones((len(xs),))
    m = len(xs)
    n = len(knots) - p - 1  # number of control points in approximation
    assert(m == len(ys))
    
    # Create matrix A
    A = np.empty((m,n))
    for row in range(m):
        for col in range(n):
            A[row, col] = np.sqrt(ws[row])*B(col, p, xs[row], knots)

    # Create vector b
    b = np.empty((len(xs),))
    for i in range(m):
        b[i] = np.sqrt(ws[i])*ys[i]

    # Compute least-squares coefficients
    Atrans = A.T
    c = np.linalg.inv(Atrans.dot(A)).dot(Atrans).dot(b)
    
    return c

def render_tensor_product_surface_alg_221(p1, p2, knots1, knots2, control_points, xs, ys):
    """ Evaluate at xs,ys """
    xs = np.array(xs)
    ys = np.array(ys)
    res = np.empty((xs.shape[0], ys.shape[0]))

    for x_i,x in enumerate(xs):
        for y_i, y in enumerate(ys):
            # Find the knot indices
            mv = get_mu(x, knots1)
            mu = get_mu(y, knots2)
            
            # Extract the correct part of the control point matrix
            C = control_points[(mv-p1):(mv+1), (mu-p2):(mu+1)]
            
            # Compute the non-zero basis functions
            Bx = alg_221(p1, knots1, x)
            By = alg_221(p2, knots2, y)
            
            value = Bx.dot(C).dot(By)
            res[x_i, y_i] = value

    return res

def B_derivative(j, p, x, knots):
    """
    Evaluate the derivative of Bj,p(x)
    p must be greater than or equal to 1
    Using theorem 3.16 from the book.
    """
    if p < 1: raise RuntimeError("p must be greater than or equal to 1")
    left = special_div(B(j, p-1, x, knots), (knots[j+p]-knots[j]) )
    right = special_div(B(j+1,p-1, x, knots), (knots[j+p+1] - knots[j+1]) )
    return (left - right)*p

def general_spline_interpolation(xs, ys, p, knots=None):
    """
    NOTE: SLOW SINCE IT USES B()
    xs,ys:  interpolation points
    p:      degree
    knots:  If None, use p+1-regular from xs[0] to slightly past x[1]
    
    returns cs, knots
    """
    
    # number of interpolation points (and also control points)
    m = len(xs)
    assert(len(ys) == m)

    # use p+1-regular knot vector with ends equal to first sample and slightly
    # past last sample
    if knots == None:
        knots = uniform_regular_knot_vector(m, p, t0=xs[0], t1=xs[-1]+0.001)

    # create matrix A
    A = np.zeros((m,m))
    for row in range(m):
        for col in range(m):
            A[row, col] = B(col, p, xs[row], knots)
    
    # compute control points
    cs = np.linalg.inv(A).dot(np.array(ys))
    return cs, knots

def alpha(j, p, i, coarse_knots, refined_knots):
    """ Compute B-splines using recursive definition. """        
    x = refined_knots[i+p]
    knots = coarse_knots
    if p == 0:
        if knots[j] <= x < knots[j+1]:
            return 1.0
        else:
            return 0.0
    else:
        left = special_div((x-knots[j])*alpha(j,p-1,i,coarse_knots, refined_knots), knots[j+p]-knots[j])
        right = special_div((knots[j+1+p]-x)*alpha(j+1,p-1,i,coarse_knots,refined_knots), knots[j+1+p]-knots[j+1])
        return left + right

def oslo_alg1(p, coarse_knots, refined_knots):
    """
    Compute knot insertion matrix using Oslo-algorithm 1.
    Returns the m x n knot insertion matrix.
    """
    n = len(coarse_knots) - p - 1
    m = len(refined_knots) - p - 1
    A = np.zeros((m, n,))
    for i in range(m):
        mu = get_mu(refined_knots[i], coarse_knots)
        if p == 0:
            res = 1
        else:
            res = R_mat(1, mu, coarse_knots, refined_knots[i+1])
            for k in range(2, p+1):
                res = res.dot(R_mat(k, mu, coarse_knots, refined_knots[i+k]))
        A[i, mu-p:(mu+1)] = res
    return A

def oslo_alg2(p, cs, coarse_knots, refined_knots):
    """
    Compute control points relative to refined knot vector
    using Oslo-algorithm 2
    Returns length-m control point vector.
    """
    n = len(coarse_knots) - p - 1
    m = len(refined_knots) - p - 1
    b = np.zeros((m,))
    for i in range(m):
        mu = get_mu(refined_knots[i], coarse_knots)
        if p == 0:
            res = cs[mu]
        else:
            _cs = cs[(mu-p):(mu+1)]
            res = R_mat(1, mu, coarse_knots, refined_knots[i+1])
            for k in range(2, p+1):
                res = res.dot(R_mat(k, mu, coarse_knots, refined_knots[i+k]))
            res = res.dot(_cs)
        b[i] = res
    return b
    
def test_render_tensor_product_surface_alg_221_1():
    """
    Compare using algorithm 2.21 with using the recurrence 
    relation.
    """
    
    n1 = 10 # x 
    n2 = 10 # y
    p1 = 3
    p2 = 3
    cs = np.random.uniform(low=0.0, high = 10.0, size=(n1, n2))
    knots1 = uniform_regular_knot_vector(n1, p1)
    knots2 = uniform_regular_knot_vector(n2, p2)
    
    # Where to evaluate
    xs = np.linspace(knots1[0], knots1[-1]-0.01, 10)
    ys = np.linspace(knots2[0], knots2[-1]-0.01, 10)
    for x in xs:
        for y in ys:
            f1 = render_tensor_product_surface_alg_221(p1, p2, knots1, knots2, cs, [x], [y])
            f2 = render_tensor_prod_spline(p1, p2, knots1, knots2, cs, [x], [y])
            assert(float_is_zero(f1-f2))
    

def test_alg_220_1():
    """
    Verify that alg 2.20 agrees with the direct implementation
    of algorithm 2.20
    """
    
    p = 3
    n = 10
    knots = uniform_regular_knot_vector(n, p)
    cs = range(10)
    for t in np.linspace(0.0, 0.99, 100):
        # Evaluate with recurrence relation
        x1 = render_spline(p, knots, cs, [t])
        # Evaluate with alg 2.20
        x2 = alg_220(p, knots, cs, t)
        assert(float_is_zero(x1-x2))

def test_get_mu_1():
    """
    Verify that the knot interval search works
    on some simple cases
    """
    assert(get_mu(0.1, [0.0, 1.0]) == 0)
    assert(get_mu(1.0, [0.0, 1.0, 2.0]) == 1)
    assert(get_mu(1.0, [0.0, 0.0, 0.0, 0.9, 1.1, 1.1, 1.1]) == 3)

def test_alg_221_1():
    """
    Verify that algorithm 2.21 produces the same
    results as the direct implementation of the
    recurrence relation.
    """
    n = 5
    p = 3
    knots = uniform_regular_knot_vector(n, p)
    
    ts = np.linspace(knots[0], knots[-1]-0.01, 10)
    
    for t in ts:
        # Compute all B-splines at once with alg 2.21
        B1 = alg_221(p, knots, t)
        mu = get_mu(t, knots)
        
        # Compute one by one with recurrence relation
        for j in range(p+1):
            B2 = B(j+mu-p, p, t, knots)
            assert( float_is_zero(B2 - B1[j]) )
        

def test_bspline_derivative():
    """
    Compare B-spline derivatives with numerical differentiation
    of the B-splines.
    """
    import matplotlib.pyplot as plt
    p = 2
    knots = [0.0, 0.0, 0.0, 1.0, 1.4, 1.9, 2.1, 2.1, 2.1]
    ts = np.linspace(knots[0], knots[-1]-0.001, 10000)
    for j in range(0, len(knots)-p-1):
        plt.figure(j)
        xs = map(lambda x: B(j, p, x, knots), ts)
        xs_der = map(lambda x: B_derivative(j, p, x, knots), ts)
        plt.subplot(2,1,1)
        plt.plot(ts, xs_der, label='Analytical')
        
        xs_num_der = np.diff(xs)/np.diff(ts)
        plt.plot(ts[:-1], xs_num_der, label='Numerical')
        plt.title('Basis function derivative %d' % j)
        plt.subplot(2,1,2)
        error = xs_der[:-1]-xs_num_der
        plt.plot(error)
        plt.title('Error')
        
    plt.show()
        
        
   
    
if __name__ == '__main__':
    test_get_mu_1()
    test_alg_220_1()
    test_alg_221_1()
    test_render_tensor_product_surface_alg_221_1()
    test_bspline_derivative()
