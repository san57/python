import copy

import numpy as np

from pycif.utils.check import info
from .aux import descentdir


def m1qn3(self, f, g, chi, **kwargs):
    # x, f, g, auxil, dxmin, df1, epsg, impres, io, mode, niter,
    #    nsim, m, jmin, jmax):
    """The optimizer itself

    Args:

       logical inmemo
       integer n,impres,io,mode,niter,nsim,m,jmin,jmax
       real x(n),f,g(n),dxmin,df1,epsg,d(n),gg(n),diag(n),aux(n),
      /     alpha(m),ybar(n,10),sbar(n,10)
       external simul,prosca,ctonb,ctcab

           variables locales

       logical sscale,cold,warm
       integer i,itmax,moderl,isim,jcour,indic
       real r1,t,tmin,tmax,gnorm,eps1,ff,preco,precos,ys,den,dk,dk1
       double precision ps,ps2,hp0

    Notes: in Version 2.0d:  July 1996:
        Wolfe's parameter (rm2 ) has been set to 0.99
        to accelerate acceptance of the step by the line-search routine
    """

    # M1QN3 Options:
    mode = self.mode
    m = self.m
    niter = self.niter
    nsim = self.nsim
    jmin = self.jmin
    jmax = self.jmax
    epsg = self.epsg
    dxmin = self.dxmin
    df1 = self.df1

    rm1 = getattr(self, "rm1", 1e-4)
    rm2 = getattr(self, "rm2", 0.99)
    rmin = getattr(self, "rmin", 1e-20)

    nverbose = kwargs.get("verbose", 2)
    logfile = kwargs.get("logfile", None)

    # M1QN3 main inputs:
    x = chi

    n = len(x)

    # Initialization
    d = np.zeros(n)
    gg = np.zeros(n)
    diag = np.ones(n)
    aux = np.zeros(n)
    alpha = np.zeros(m)
    ybar = np.zeros((10, n))
    sbar = np.zeros((10, n))

    sscale = 1
    if mode - int(mode / 2.0) * 2.0 == 0:
        sscale = 0

    warm = 0
    if mode // 2 == 1:
        warm = 1

    cold = not warm

    itmax = niter
    niter = 0
    isim = 1
    eps1 = 1.0

    ps = np.dot(g, g)
    dk = ps
    gnorm = ps
    gnorm = np.sqrt(gnorm)
    if nverbose >= 1:
        info("     f         = " + str(f))
        info("     norm of g = " + str(gnorm))
    if gnorm < rmin:
        mode = 2
        if nverbose >= 1:
            info(" >>> m1qn3a: initial gradient is too small")
        nsim = isim
        epsg = eps1
        return x, f, g, niter, nsim, epsg, mode
    #
    #     --- initializing descent direction
    #
    if cold:
        jmin = 0
        jmax = -1
    jcour = jmax
    #
    #     --- mise a l'echelle de la premiere direction de descente
    #
    if cold:
        #
        #         --- use Fletcher's scaling and initialize diag to 1.
        #
        precos = 2.0 * df1 / (gnorm * gnorm)
        d = -g * precos
        if nverbose >= 5:
            info(" m1qn3a: descent direction -g: precon = " + str(precos))
    else:
        #
        #         --- use the matrix stored in [diag and] the (y,s) pairs
        #
        if sscale:
            ps = np.dot(ybar[jcour, :], ybar[jcour, :])
            precos = 1.0 / ps

        ctonb = None
        ctcab = None

        d = descentdir(
            n, sscale, m, -g, jmin, jmax, precos, diag, alpha, ybar, sbar
        )

    #
    #     --- initialisation pour mlis0
    #
    tmax = 1.0e20
    hp0 = np.dot(d, g)
    if hp0 >= 0.0:
        mode = 7
        if nverbose >= 1:
            info(" >>> m1qn3 (iteration " + str(niter) + "):")
            info(
                "     the search direction d is not a descent direction: (g,"
                "d) = " + str(hp0)
            )
            nsim = isim
            epsg = eps1
            return x, f, g, niter, nsim, epsg, mode
    #
    #     --- compute the angle (-g,d)
    #
    if warm and nverbose >= 5:
        ps = np.sqrt(np.dot(g, g))
        ps2 = np.sqrt(np.dot(d, d))
        ps = hp0 / ps / ps2
        ps = min(-ps, 1.0)
        r1 = np.degrees(np.arccos(ps))
        info(
            " m1qn3: descent direction d: angle(-g,d) = "
            + str(r1)
            + "degrees",
        )
    #
    # ---- Debut de l'iteration. on cherche x(k+1) de la forme x(k) + t*d,
    #     avec t > 0. On connait d.
    #
    #         Debut de la boucle: etiquette 100,
    #         Sortie de la boucle: goto 1000.
    #
    while 1:
        niter += 1
        #   if verbose < 0:
        #     if niter - int(niter/(-verbose))*(-verbose) == 0 :
        #       indic=1
        #       f, g = simul.simul (x,niter+999,auxil,io)
        #       continue
        if nverbose >= 3:
            info(
                " m1qn3: iter "
                + str(niter)
                + " simul "
                + str(isim)
                + ", f="
                + str(f)
                + ", h'(0)="
                + str(hp0)
            )
        gg = copy.copy(g)
        ff = f
        #
        #     --- recherche lineaire et nouveau point x(k+1)
        #
        if nverbose >= 5:
            info(" m1qn3: line search")
        #
        #         --- calcul de tmin
        #
        tmin = 0.0
        for i in range(n):
            tmin = max(tmin, abs(d[i]))
        tmin = dxmin / tmin
        t = 1.0
        r1 = hp0

        x, f, g, isim, moderl = self.mlis0(
            x,
            f,
            r1,
            t,
            tmin,
            tmax,
            d,
            g,
            rm2,
            rm1,
            nverbose,
            logfile,
            isim,
            nsim,
            niter,
            **kwargs
        )
        #
        #         --- mlis0 renvoie les nouvelles valeurs de x, f et g
        #
        if moderl != 0:
            if moderl < 0:
                #
                #             --- calcul impossible
                #                 t, g: ou les calculs sont impossibles
                #                 x, f: ceux du t_gauche (donc f <= ff)
                #
                mode = moderl
            elif moderl == 1:
                #
                #             --- descente bloquee sur tmax
                #                 [sortie rare (!!) d'apres le code de mlis0]
                #
                mode = 3
                if nverbose >= 1:
                    info(
                        " >>> m1qn3 (iteration "
                        "): line search blocked on tmax: decrease "
                        "the scaling".format(niter)
                    )
            elif moderl == 4:
                #
                #             --- nsim atteint
                #                 x, f: ceux du t_gauche (donc f <= ff)
                #
                mode = 5
            elif moderl == 5:
                #
                #             --- arret demande par l'utilisateur (indic = 0)
                #                 x, f: ceux en sortie du simulateur
                #
                mode = 0
            elif moderl == 6:
                #
                #             --- arret sur dxmin ou appel incoherent
                #                 x, f: ceux du t_gauche (donc f <= ff)
                #
                mode = 6
            nsim = isim
            epsg = eps1
            return x, f, g, niter, nsim, epsg, mode
        #
        # NOTE: stopping tests are now done after having updated the matrix, so
        # that update information can be stored in case of a later warm restart
        #
        #     --- mise a jour de la matrice
        #
        if m > 0:
            #
            #         --- mise a jour des pointeurs
            #
            jmax += 1
            if jmax >= m:
                jmax = jmax - m

            if (cold and niter > m) or (warm and jmin == jmax):
                jmin += 1
                if jmin >= m:
                    jmin -= m
            jcour = jmax
            #
            #         --- y, s et (y,s)
            #
            sbar[jcour, :] = t * d
            ybar[jcour, :] = g - gg
            if nverbose >= 5:
                ps = np.dot(sbar[jcour, :], sbar[jcour, :])
                dk1 = np.sqrt(ps)
                if niter > 1:
                    info(
                        " m1qn3: convergence rate, s(k)/s(k-1) = "
                        + str(dk1 / dk)
                    )
                dk = dk1
            ps = np.dot(ybar[jcour, :], sbar[jcour, :])
            ys = ps
            if ys <= 0.0:
                mode = 7
                if nverbose >= 1:
                    info(
                        " >>> m1qn3 (iteration "
                        + str(niter)
                        + "): the scalar product (y,s) = "
                        + str(ys)
                        + "is not positive"
                    )
                nsim = isim
                epsg = eps1
                return x, f, g, niter, nsim, epsg, mode
            #
            #         --- ybar et sbar
            #
            r1 = np.sqrt(1.0 / ys)
            sbar[jcour, :] = r1 * sbar[jcour, :]
            ybar[jcour, :] = r1 * ybar[jcour, :]
            #
            #         --- compute the scalar or diagonal preconditioner
            #
            if nverbose >= 5:
                info(" m1qn3: matrix update:")
            #
            #             --- Here is the Oren-Spedicato factor, for scalar
            # scaling
            #
            if sscale:
                ps = np.dot(ybar[jcour, :], ybar[jcour, :])
                precos = 1.0 / ps

                if nverbose >= 5:
                    info("Oren-Spedicato factor = " + str(precos))
            #
            #             --- Scale the diagonal to Rayleigh's ellipsoid.
            #                 Initially (niter.eq.1) and for a cold start,
            # this is
            #                 equivalent to an Oren-Spedicato scaling of the
            #                 identity matrix.
            #
            else:
                aux = ybar[jcour, :]
                ps = 0.0
                for i in range(n):
                    ps += diag[i] * aux[i] * aux[i]
                r1 = 1.0 / ps
                if nverbose >= 5:
                    info("     fitting the ellipsoid: factor = " + str(r1))
                diag = diag * r1
                #
                #             --- update the diagonal
                #                 (gg is used as an auxiliary vector)
                #
                gg = sbar[jcour, :]
                ps = 0.0
                for i in range(n):
                    ps += gg[i] * gg[i] / diag[i]

                den = ps
                iprint = 0
                for i in range(n):
                    diag[i] = 1.0 / (
                        1.0 / diag[i]
                        + aux[i] * aux[i]
                        - (gg[i] / diag[i]) * (gg[i] / diag[i]) / den
                    )
                    if diag[i] <= 0.0:
                        diag[i] = rmin
                        iprint = i
                if iprint != 0 and nverbose >= 5:
                    info(
                        " >>> m1qn3-WARNING: diagonal element "
                        + str(iprint)
                        + " is negative ("
                        + str(diag[iprint])
                        + "), reset to "
                        + str(rmin)
                    )

                if nverbose >= 5:
                    ps = 0.0
                    for i in range(n):
                        ps += diag[i]
                    ps = ps / n
                    preco = ps

                    ps2 = 0.0
                    for i in range(n):
                        ps2 += (diag[i] - ps) * (diag[i] - ps)

                    ps2 = np.sqrt(ps2 / n)
                    info(
                        "     updated diagonal: average value = "
                        + str(preco)
                        + ", sqrt(variance) = "
                        + str(ps2)
                    )
        #
        #     --- tests d'arret
        #
        ps = np.dot(g, g)
        eps1 = ps
        eps1 = np.sqrt(eps1) / gnorm

        # Some verbose
        toprint = """
        M1QN3:
        Iteration number: {}
        Simulation number: {} over {}
        Cost function: {}
        Gradient norm ratio: {}
        """.format(
            niter, isim, nsim, f, eps1
        )
        info(toprint)

        if nverbose >= 5:
            info(" m1qn3: stopping criterion on g: " + str(eps1))
        if eps1 < epsg:
            mode = 1
            nsim = isim
            epsg = eps1
            return x, f, g, niter, nsim, epsg, mode
        if niter == itmax:
            mode = 4
            if nverbose >= 1:
                info(
                    " >>> m1qn3 (iteration "
                    + str(niter)
                    + "): maximal number of iterations"
                )
            nsim = isim
            epsg = eps1
            return x, f, g, niter, nsim, epsg, mode
        if isim > nsim:
            mode = 5
            if nverbose >= 1:
                info(
                    " >>> m1qn3 (iteration "
                    + str(niter)
                    + "): "
                    + str(isim)
                    + " simulations (maximal number reached)"
                )
            nsim = isim
            epsg = eps1
            return x, f, g, niter, nsim, epsg, mode
        #
        #     --- calcul de la nouvelle direction de descente d = - H.g
        #
        if m == 0:
            preco = 2.0 * (ff - f) / ((eps1 * gnorm) * (eps1 * gnorm))
            d[:] = -g[:] * preco
        else:
            info("Computing new descent direction")
            d = descentdir(
                n, sscale, m, -g, jmin, jmax, precos, diag, alpha, ybar, sbar
            )
        #
        #         --- test: la direction d est-elle de descente ?
        #             hp0 sera utilise par mlis0
        #
        hp0 = np.dot(d, g)
        if hp0 >= 0.0:
            mode = 7
            if nverbose >= 1:
                info(" >>> m1qn3 (iteration " + str(niter) + "):")
                info(
                    "     the search direction d is not a descent direction: "
                    "(g,d) = " + str(hp0)
                )
            nsim = isim
            epsg = eps1
            return x, f, g, niter, nsim, epsg, mode

        if nverbose >= 5:
            ps = np.dot(g, g)
            ps = np.sqrt(ps)
            ps2 = np.dot(d, d)
            ps2 = np.sqrt(ps2)
            ps = hp0 / ps / ps2
            ps = min(-ps, 1.0)
            r1 = np.degrees(np.arccos(ps))
            info(
                " m1qn3: descent direction d: angle(-g,d) = "
                + str(r1)
                + "degrees"
            )
