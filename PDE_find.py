import numpy as np


def FiniteDiff(u, dx):

    n = u.size
    ux = np.zeros(n, dtype=np.complex64)

    for i in range(1, n - 1):
        ux[i] = (u[i + 1] - u[i - 1]) / (2 * dx)

    ux[0] = (-3.0 / 2 * u[0] + 2 * u[1] - u[2] / 2) / dx
    ux[n - 1] = (3.0 / 2 * u[n - 1] - 2 * u[n - 2] + u[n - 3] / 2) / dx
    return ux


def Diff(u, dxt, name):
    """
    Here dx is a scalar, name is a str indicating what it is
    """

    n, m = u.shape
    uxt = np.zeros((n, m), dtype=np.complex64)

    if name == 'x':
        for i in range(m):
            uxt[:, i] = FiniteDiff(u[:, i], dxt)

    elif name == 't':
        for i in range(n):
            uxt[i, :] = FiniteDiff(u[i, :], dxt)

    else:
        NotImplementedError()

    return uxt


def Train(R, Ut, lam, d_tol, maxit=10, STR_iters=10, l0_penalty=1, normalize=2, split=0.8,
          print_best_tol=False, sparse='STR'):
    """
    This function trains a predictor using STRidge.
    It runs over different values of tolerance and trains predictors on a training set, then evaluates them
    using a loss function on a holdout set.
    """

    # Split data into 80% training and 20% test, then search for the best tolderance.
    # np.random.seed(0)  # for consistancy
    n, _ = R.shape
    #train = np.random.choice(n, int(n * split), replace=False)
    #test = [i for i in np.arange(n) if i not in train]
    TrainR = R#[train, :]
    TestR = R#[test, :]
    TrainY = Ut#[train, :]
    TestY = Ut#[test, :]
    D = TrainR.shape[1]

    # Set up the initial tolerance and l0 penalty
    d_tol = float(d_tol)
    tol = d_tol
    if l0_penalty == None: l0_penalty = 0.001 * np.linalg.cond(R)

    # Get the standard least squares estimator
    w = np.zeros((D, 1))

    # check
    # print(np.nan in TrainR)
    # print(np.inf in TrainR)

    def AIC(w, err):
        k = 0
        for item in w:
            if item != 0:
                k += 1

        return 2*k+2*np.log(err)

    w_best = np.linalg.lstsq(TrainR, TrainY)[0]
    data_err_best = np.linalg.norm(TestY - TestR.dot(w_best), 2)
    err_best = np.linalg.norm(TestY - TestR.dot(w_best), 2) + l0_penalty * np.count_nonzero(w_best)
    aic_best = AIC(w_best[:, 0], data_err_best)
    tol_best = 0

    if sparse == 'STR':
        # Now increase tolerance until test performance decreases
        for iter in range(maxit):

            # Get a set of coefficients and error
            w = STRidge(R, Ut, lam, STR_iters, tol, normalize=normalize)
            err = np.linalg.norm(TestY - TestR.dot(w), 2) + l0_penalty * np.count_nonzero(w)
            data_err = np.linalg.norm(TestY - TestR.dot(w), 2)

            # Has the accuracy improved?
            aic = AIC(w[:, 0], data_err)
            if aic <= aic_best:
                aic_best = aic
                err_best = err
                w_best = w
                data_err_best = data_err
                tol_best = tol
                tol = tol + d_tol
            else:
                tol = max([0, tol - 2 * d_tol])
                d_tol = 2 * d_tol / (maxit - iter)
                tol = tol + d_tol

        if print_best_tol: print("Optimal tolerance:", tol_best)

    elif sparse == 'Lasso':
        w = Lasso(R, Ut, lam, w=np.array([0]), maxit=maxit*10, normalize=normalize)
        err = np.linalg.norm(Ut - R.dot(w), 2) + l0_penalty * np.count_nonzero(w)
        data_err = np.linalg.norm(Ut - R.dot(w), 2)

        if err <= err_best:
            err_best = err
            w_best = w
            data_err_best = data_err

    return w_best, err_best, data_err_best, aic_best


def STRidge(X0, y, lam, maxit, tol, normalize=0, print_results=False):
    """
    Sequential Threshold Ridge Regression algorithm for finding (hopefully) sparse
    approximation to X^{-1}y.  The idea is that this may do better with correlated observables.
    This assumes y is only one column
    """

    n, d = X0.shape
    X = np.zeros((n, d), dtype=np.complex64)
    # First normalize data
    if normalize != 0:
        Mreg = np.zeros((d, 1))
        for i in range(0, d):
            Mreg[i] = 1.0 / (np.linalg.norm(X0[:, i], normalize))
            X[:, i] = Mreg[i] * X0[:, i]
    else:
        X = X0

    # Get the standard ridge estimate
    if lam != 0:
        w = np.linalg.lstsq(X.T.dot(X) + lam * np.eye(d), X.T.dot(y))[0]
    else:
        w = np.linalg.lstsq(X, y)[0]
    num_relevant = d
    biginds = np.where(abs(w) > tol)[0]

    # Threshold and continue
    for j in range(maxit):

        # Figure out which items to cut out
        smallinds = np.where(abs(w) < tol)[0]
        new_biginds = [i for i in range(d) if i not in smallinds]

        # If nothing changes then stop
        if num_relevant == len(new_biginds):
            break
        else:
            num_relevant = len(new_biginds)

        # Also make sure we didn't just lose all the coefficients
        if len(new_biginds) == 0:
            if j == 0:
                # if print_results: print "Tolerance too high - all coefficients set below tolerance"
                return w
            else:
                break
        biginds = new_biginds

        # Otherwise get a new guess
        w[smallinds] = 0
        if lam != 0:
            w[biginds] = \
            np.linalg.lstsq(X[:, biginds].T.dot(X[:, biginds]) + lam * np.eye(len(biginds)), X[:, biginds].T.dot(y))[0]
        else:
            w[biginds] = np.linalg.lstsq(X[:, biginds], y)[0]

    # Now that we have the sparsity pattern, use standard least squares to get w
    if biginds != []: w[biginds] = np.linalg.lstsq(X[:, biginds], y)[0]

    if normalize != 0:
        return np.multiply(Mreg, w)
    else:
        return w


def Lasso(X0, Y, lam, w=np.array([0]), maxit=100, normalize=2):
    """
    Uses accelerated proximal gradient (FISTA) to solve Lasso
    argmin (1/2)*||Xw-Y||_2^2 + lam||w||_1
    """

    # Obtain size of X
    n, d = X0.shape
    X = np.zeros((n, d), dtype=np.complex64)
    Y = Y.reshape(n, 1)

    # Create w if none is given
    if w.size != d:
        w = np.zeros((d, 1), dtype=np.complex64)
    w_old = np.zeros((d, 1), dtype=np.complex64)

    # First normalize data
    if normalize != 0:
        Mreg = np.zeros((d, 1))
        for i in range(0, d):
            Mreg[i] = 1.0 / (np.linalg.norm(X0[:, i], normalize))
            X[:, i] = Mreg[i] * X0[:, i]
    else:
        X = X0

    # Lipschitz constant of gradient of smooth part of loss function
    L = np.linalg.norm(X.T.dot(X), 2)

    # Now loop until converged or max iterations
    for iters in range(0, maxit):

        # Update w
        z = w + iters / float(iters + 1) * (w - w_old)
        w_old = w
        z = z - X.T.dot(X.dot(z) - Y) / L
        for j in range(d):
            w[j] = np.multiply(np.sign(z[j]), np.max([abs(z[j]) - lam / L, 0]))

        # Could put in some sort of break condition based on convergence here.

    # Now that we have the sparsity pattern, used least squares.
    biginds = np.where(w != 0)[0]
    if biginds != []: w[biginds] = np.linalg.lstsq(X[:, biginds], Y)[0]

    # Finally, reverse the regularization so as to be able to use with raw data
    if normalize != 0:
        return np.multiply(Mreg, w)
    else:
        return w


def build_system(u, dt, dx, D=4, C=1, time_diff='FD', space_diff='FD', width_x=None,
                      width_t=None, deg_x=5, deg_t=None):
    """
    Constructs a large linear system to use in later regression for finding PDE.
    This function works when we are not subsampling the data or adding in any forcing.
    Input:
        Required:
            u = data to be fit to a pde
            dt = temporal grid spacing
            dx = spatial grid spacing
        Optional:
            D = max derivative to include in rhs (default = 3)
            C = degree of polynomials to the derivative terms
            time_diff = method for taking time derivative
                        options = 'poly', 'FD', 'FDconv','TV'
                        'poly' (default) = interpolation with polynomial
                        'FD' = standard finite differences
                        'FDconv' = finite differences with convolutional smoothing
                                   before and after along x-axis at each timestep
                        'Tik' = Tikhonov (takes very long time)
            space_diff = same as time_diff with added option, 'Fourier' = differentiation via FFT
            lam_t = penalization for L2 norm of second time derivative
                    only applies if time_diff = 'TV'
                    default = 1.0/(number of timesteps)
            lam_x = penalization for L2 norm of (n+1)st spatial derivative
                    default = 1.0/(number of gridpoints)
            width_x = number of points to use in polynomial interpolation for x derivatives
                      or width of convolutional smoother in x direction if using FDconv
            width_t = number of points to use in polynomial interpolation for t derivatives
            deg_x = degree of polynomial to differentiate x
            deg_t = degree of polynomial to differentiate t
            sigma = standard deviation of gaussian smoother
                    only applies if time_diff = 'FDconv'
                    default = 2
    Output:
        ut = column vector of length u.size
        R = matrix with ((D+1)*(P+1)) of column, each as large as ut
        rhs_description = description of what each column in R is
    """

    n, m = u.shape

    if width_x == None: width_x = n / 10
    if width_t == None: width_t = m / 10
    if deg_t == None: deg_t = deg_x

    # If we're using polynomials to take derviatives, then we toss the data around the edges.
    if time_diff == 'poly':
        m2 = m - 2 * width_t
        offset_t = width_t
    else:
        m2 = m
        offset_t = 0
    if space_diff == 'poly':
        n2 = n - 2 * width_x
        offset_x = width_x
    else:
        n2 = n
        offset_x = 0

    ########################
    # First take the time derivaitve for the left hand side of the equation
    ########################
    ut = np.zeros((n2, m2), dtype=np.complex64)

    if time_diff == 'poly':
        T = np.linspace(0, (m - 1) * dt, m)
        for i in range(n2):
            ut[i, :] = PolyDiff(u[i + offset_x, :], T, diff=1, width=width_t, deg=deg_t)[:, 0]

    else:
        for i in range(n2):
            ut[i, :] = FiniteDiff(u[i + offset_x, :], dt, 1)

    ut = np.reshape(ut, (n2 * m2, 1), order='F')

    ########################
    # Now form the rhs one column at a time, and record what each one is
    ########################

    u2 = u[offset_x:n - offset_x, offset_t:m - offset_t]
    Theta = np.zeros((n2 * m2, (D + 1) * C), dtype=np.complex64)
    ux = np.zeros((n2, m2), dtype=np.complex64)
    rhs_description = ['' for i in range((D + 1) * C)]

    if space_diff == 'poly':
        Du = {}
        for i in range(m2):
            Du[i] = PolyDiff(u[:, i + offset_t], np.linspace(0, (n - 1) * dx, n), diff=D, width=width_x, deg=deg_x)
    if space_diff == 'Fourier': ik = 1j * np.fft.fftfreq(n) * n

    for d in range(D + 1):
        # compute derivatives of d degree
        if d > 0:
            for i in range(m2):
                if space_diff == 'FD':
                    ux[:, i] = FiniteDiff(u[:, i + offset_t], dx, d)
                elif space_diff == 'poly':
                    ux[:, i] = Du[i][:, d - 1]
        else:
            ux = np.array(u2, dtype=np.complex64)
        # if d == 1: print(ux)

        # compute polynomials of all terms, c used as c+1
        for c in range(C):
            Theta[:, d * C + c] = np.reshape(np.power(ux, c+1), (n2 * m2), order='F')
            # print('d:{}, c:{}, mean:{}'.format(d, c, np.mean(Theta[:, d * C + c])))

            if d > 0:
                rhs_description[d * C + c] = rhs_description[d * C + c] + \
                                             'u_{' + ''.join(['x' for _ in range(d)]) + '}'
            else:
                rhs_description[d * C + c] = rhs_description[d * C + c] + 'u'

            if c > 0:
                rhs_description[d * C + c] = rhs_description[d * C + c] + '^' + str(c+1)

    # print(rhs_description)
    features, rhs = create_cross_features(Theta, rhs_description)
    features = np.concatenate((Theta, features), 1)
    rhs = np.concatenate((rhs_description, rhs), 0)

    return ut, features, rhs


def create_cross_features(features, des):
    from itertools import combinations
    comb_f = list(combinations(features.transpose((1,0)), 2))
    comb_d = list(combinations(des, 2))
    cross_features = []
    cross_des = []
    for item in comb_f:
        cross_features.append(np.array(item[0]) * np.array(item[1]))
    for item in comb_d:
        cross_des.append(item[0] + '*' + item[1])
    cross_features = np.array(cross_features)

    return cross_features.transpose((1,0)), cross_des


def eq_pde(w, rhs_description, ut='u_t'):
    pde = ut + ' = '
    first = True
    for i in range(len(w)):
        if abs(w[i]) >= 1e-4:
            if not first:
                pde = pde + ' + '
            pde = pde + "%05f" % w[i].real + rhs_description[i] + "\n   "
            first = False
    return pde