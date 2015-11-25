import bsplines

def compute_num_basis_funs(p, knots):
    # compute number of basis functions from the length
    # of the knot vector using relation Nk = n + p + 1.
    num_basis_funs = len(knots)-p-1
    return num_basis_funs
    
def b(j, p, x, knots, weights):
    """
    Evaluate rational B-spline basis function no. j,
    of degree p with knot vector knots, at the point x.
    """
    num_basis_funs = compute_num_basis_funs(p, knots)

    assert(len(weights) == num_basis_funs)
    den = 0.0
    for i in range(num_basis_funs):
        den += bsplines.B(i, p, x, knots)*weights[i]
    num = bsplines.B(j, p, x, knots)*weights[j]
    return bsplines.special_div(num, den)

def b_derivative(j, p, x, knots, weights):
    """
    Derivative of rational B-spline basis functions.
    """
    num_basis_funs = compute_num_basis_funs(p, knots)
    num_basis_funs = len(weights)
    assert(len(weights) == num_basis_funs)
    
    weighted_sum = 0.0
    for i in range(num_basis_funs):
        weighted_sum += bsplines.B(i, p, x, knots)*weights[i]
    
    left_num = weighted_sum*bsplines.B_derivative(j, p, x, knots)*weights[j]
    
    right_num = 0.0
    for i in range(num_basis_funs):
        right_num += bsplines.B_derivative(i, p, x, knots)*weights[i]
    right_num *= bsplines.B(j, p, x, knots)*weights[j]

    den = weighted_sum*weighted_sum
    return bsplines.special_div(left_num-right_num, den)

def render_spline(p, knots, control_points, weights, ts):
    """
    Compute points on a spline function using the straightforward
    implementation of the recurrence relation for the NURBS basis
    function.
    """
    assert len(control_points) == len(weights)
    ys = []
    for t in ts:
        y = 0.0 
        for j in range(0, len(control_points)):
            y += b(j, p, t, knots, weights)*control_points[j]
        ys.append(y)
    return ys  

def render_spline_derivative(p, knots, control_points, weights, ts):
    """
    Compute derivative point on a spline function using the straightforward
    implementation of the recurrence relation for the NURBS basis
    function.
    """
    assert len(control_points) == len(weights)
    ys = []
    for t in ts:
        y = 0.0 
        for j in range(0, len(control_points)):
            y += b_derivative(j, p, t, knots, weights)*control_points[j]
        ys.append(y)
    return ys  

    
def test1():
    import matplotlib.pyplot as plt
    import numpy as np
    
    xs = [1, 1, 0, -1, -1, -1, 0, 1]
    ys = [0, 1, 1, 1, 0, -1, -1, -1]
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    degree = 3
    knots = [0.0, 0.0, 0.0, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.0, 1.0, 1.01]
    
    # Render curve
    ts = np.linspace(0.0, 1.0, 100)
    xs_vis = []
    ys_vis = []
    for t in ts:
        x = 0.0
        y = 0.0
        for i in range(len(xs)):
            temp = b(i, degree, t, knots, weights)
            x += temp*xs[i]
            y += temp*ys[i]
        xs_vis.append(x)
        ys_vis.append(y)
    
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')
    
    ax.plot(xs, ys, label='control grid', marker='o')
    ax.plot(xs_vis, ys_vis, label='model')

    # Render normals
    length = 0.5
    ts_normal = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.9, 0.9999]
    for t in ts_normal:
        x = 0.0
        y = 0.0
        dx = 0.0
        dy = 0.0
        for i in range(len(xs)):
            temp = b(i, degree, t, knots, weights)
            x += temp*xs[i]
            y += temp*ys[i]
            temp_derivative = b_derivative(i, degree, t, knots, weights)
            dx += temp_derivative*xs[i]
            dy += temp_derivative*ys[i]
        tangent = np.array([dx,dy,0.0])
        z_axis = np.array([0.0, 0.0, 1.0])
        normal = np.cross(tangent,z_axis)
        normal = normal / np.linalg.norm(normal)
        
        print normal
        
        ax.plot([x, x+normal[0]*length], [y, y+normal[1]*length])

    plt.legend()
    plt.show()

def test2():
    """ Verify NURBS unit circle representation. """
    import matplotlib.pyplot as plt
    import numpy as np
    fig = plt.figure()
    fig.add_subplot(111, aspect='equal')

    w = np.sqrt(2)/2
    xs = [0.0,  -1.0, -1.0, -1.0, 0.0, 1.0, 1.0, 1.0,  0.0]
    ys = [-1.0, -1.0, 0.0,  1.0,  1.0, 1.0, 0.0, -1.0, -1.0]
    ws = [1,    w,    1,    w,    1,   w,   1,   w,    1]
    p = 2
    knots = [0,0,0,0.25,0.25,0.5,0.5,0.75,0.75,1,1,1]
    
    ts = np.linspace(0.0, 0.99, 100)
    xs_vis = []
    ys_vis = []
    xs_der_vis = []
    ys_der_vis = []
    for t in ts:
        x = 0.0
        y = 0.0
        x_der = 0.0
        y_der = 0.0
        for j in range(len(xs)):
            eval_basis = b(j, p, t, knots, ws)
            eval_basis_der = b_derivative(j, p, t, knots, ws)
            x += eval_basis*xs[j]
            y += eval_basis*ys[j]
            x_der += eval_basis_der*xs[j]
            y_der += eval_basis_der*ys[j]
        xs_vis.append(x)
        ys_vis.append(y)
        xs_der_vis.append(x_der)
        ys_der_vis.append(y_der)
    
    plt.plot(xs, ys, marker='x')
    plt.plot(xs_vis, ys_vis)
    
    # draw the curve normals
    z_ax = np.array([0.0, 0.0, 1.0])
    for i in range(len(xs_vis)):
        tangent = np.array([xs_der_vis[i],ys_der_vis[i],0.0])
        normal = np.cross(z_ax, tangent)
        normal = normal / np.linalg.norm(normal)
        l = -1
        plt.plot([xs_vis[i], xs_vis[i]+l*normal[0]], [ys_vis[i], ys_vis[i]+l*normal[1]], c='k')
    
        # Verify that each normal passes through the origin using distance formula
        # from point x0 to line passing through x1 and x2.
        x1 = np.array([xs_vis[i], ys_vis[i], 0.0])
        x2 = np.array([xs_vis[i]+normal[0], ys_vis[i]+normal[1], 0.0])
        x0 = np.array([0,0,0])
        dist = np.linalg.norm(np.cross(x0-x1, x0-x2)) / np.linalg.norm(x2-x1)
        assert dist < 1e-6
    
    # verify that each vis point is on the unit circle
    for x,y in zip(xs_vis, ys_vis):
        assert abs( x**2 + y**2 - 1) < 1e-6
    
    
    plt.show()
    
    

if __name__ == '__main__':
    test1()
    test2()        
