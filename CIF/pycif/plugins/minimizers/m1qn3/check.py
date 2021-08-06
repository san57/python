import numpy as np

from pycif.utils.check import info


def check_options(self, chi, finit, **kwargs):
    """Check the consistency of M1QN3 parameters and initializes mode and
    working zone

    Args:
        chi (vector): initial value for the vector to optimize
        finit (float): initial value of the function
        minimizer (Plugin): Options for the minimizer

    Kwargs:
        dim (int): dimension of the problem
        niter (int): number of iterations before stop
        nsim (int): number of simulations before stop
        iz (list): list of working arguments
        df1 (float): expected decrease for f
        m (int): number of updates
        dxmin (float): precision for x. Default is 1e-20
        epsg (float): relative precision of g. default is 1e-20
        impres (int): verbose level
        mode (int): mode for M1QN3
        logfile (str): file for verbose entries

    Return:
        (mode, iz): tuple with the updated model and working zone updated in
        kwargs

    """

    # Dimension of the control vector
    dim = chi.size

    # Number of simulations
    # Can be explicitly given by the user, or will be determined from 'maxiter'
    try:
        niter = self.niter
        nsim = self.nsim

    except AttributeError:
        maxiter = self.maxiter
        niter = maxiter
        nsim = 2 * maxiter

    # Default parameters
    m = getattr(self, "m", 5)
    dxmin = getattr(self, "dxmin", 1.0e-20)
    epsg = getattr(self, "epsg", 1.0e-20)
    mode = getattr(self, "mode", 0)
    df1 = getattr(self, "df1", 0.01) * finit
    iz = getattr(self, "iz", np.zeros(6, dtype=int))

    # Init parameters
    towrite = """
        M1QN3 (Version 2.0b, December 1993):
        entry point program translated in Python language
           dimension of the problem (n): {}
           number of updates: {}
           absolute precision on x (dxmin) : {}
           expected decrease for f (df1): {}
           relative precision on g (epsg): {}
           maximal number of iterations (niter): {}
           maximal number of simulations (nsim): {}
        """.format(
        dim, m, dxmin, df1, epsg, niter, nsim
    )

    info(towrite)

    # Checking for inconsistent definition of parameters
    if (
        dim <= 0.0
        or niter <= 0.0
        or nsim <= 0.0
        or dxmin <= 0.0
        or epsg <= 0.0
        or epsg > 1.0
        or mode < 0.0
        or mode > 3
        or m < 1
        or m > 10
    ):
        info("Warning: inconsistent call in M1QN3")
        raise ValueError

    # Select mode
    if mode - int(mode / 2.0) * 2.0 == 0:
        info("M1QN3: Diagonal Initial Scaling mode")
        sscale = False

    else:
        info("M1QN3: Scalar Initial Scaling mode")
        sscale = True

    # Cold start or warm restart?
    # Check iz: iz(1)=n, iz(2)=(0 if DIS, 1 if SIS),
    #           iz(3)=m, iz(4)=jmin, iz(5)=jmax
    if mode % 2 == 0:
        info("M1QN3: cold start")

    else:
        if (
            iz[1] != dim
            or iz[2] != sscale
            or iz[3] != m
            or iz[4] < 1
            or iz[5] < 1
            or iz[4] > iz[3]
            or iz[5] > iz[3]
        ):
            info("Warning: M1QN3: inconsistent iz for a warm " "restart")
            raise ValueError
        info("M1QN3: warm restart")

    iz[1] = dim
    iz[2] = 0
    if sscale:
        iz[2] = 1

    iz[3] = m

    # Update options
    self.mode = mode
    self.iz = iz
    self.m = iz[3]
    self.jmin = iz[4]
    self.jmax = iz[5]
    self.niter = niter
    self.nsim = nsim
    self.dxmin = dxmin
    self.df1 = df1
    self.epsg = epsg

    return self
