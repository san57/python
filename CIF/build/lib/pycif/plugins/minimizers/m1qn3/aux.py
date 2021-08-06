from __future__ import division

import numpy as np
from builtins import range
from builtins import str

from pycif.utils.check import info


def mlis0(
    self,
    xn,
    fn,
    fpn,
    t,
    tmin,
    tmax,
    d,
    g,
    amd,
    amf,
    imp,
    io,
    nap,
    napmax,
    niter,
    **kwargs
):
    """

       mlis0 + minuscules + commentaires
       + version amelioree (XII 88): interpolation cubique systematique
         et anti-overflows
       + declaration variables (II/89, JCG).
       + barr is also progressively decreased (12/93, CL & JChG).
         barmul is set to 5.

       ----------------------------------------------------------------

          en sortie logic =

          0          descente serieuse
          1          descente bloquee
          4          nap > napmax
          5          retour a l'utilisateur
          6          fonction et gradient pas d'accord
          < 0        contrainte implicite active

    ----

    --- arguments

       external simul,prosca
       integer n,imp,io,logic,nap,napmax
       real xn(n),fn,fpn,t,tmin,tmax,d(n),g(n),amd,amf,x(n)

    --- variables locales

       integer i,indic,indica,indicd
       real tesf,tesd,tg,fg,fpg,td,ta,fa,fpa,d2,f,fp,ffn,fd,fpd,
      1 z,test,barmin,barmul,barmax,barr,gauche,droite,
      2 taa
       double precision ps
    """

    # Simulator
    simulator = self.simulator

    n = len(xn)

    if not (
        n > 0
        and fpn < 0.0
        and t > 0.0
        and tmax > 0.0
        and amf > 0.0
        and amd > amf
        and amd < 1.0
    ):
        toprint = """
        MODE == 6!!!!
        n = {}
        fpn = {}
        t = {}
        tmax = {}
        amf = {}
        amd = {}
        """.format(
            n, fpn, t, tmax, amf, amd
        )
        info(toprint, io)

        logic = 6
        return xn, fn, g, nap, logic

    tesf = amf * fpn
    tesd = amd * fpn
    barmin = 0.01
    barmul = 5.0
    barmax = 0.3
    barr = barmin
    td = 0.0
    tg = 0.0
    fg = fn
    fpg = fpn
    ta = 0.0
    fa = fn
    fpa = fpn
    ps = np.dot(d, d)
    d2 = ps
    #
    #               elimination d'un t initial ridiculement petit
    #
    t = max(t, tmin)
    if t > tmax:
        tmin = tmax
        info("     mlis0          tmin forced to tmax", io)

    while (fn + t * fpn) >= (fn + 0.9 * t * fpn):
        t = 2.0 * t

    indica = 1
    logic = 0
    if t > tmax:
        t = tmax
        logic = 1

    if imp >= 4:
        info(
            "     mlis0   fpn="
            + str(fpn)
            + " d2="
            + str(d2)
            + " tmin="
            + str(tmin)
            + " tmax="
            + str(tmax)
            + "degrees",
            io,
        )

    #
    #     --- nouveau x
    #
    x = xn + t * d

    #
    # --- boucle
    #
    while 1:

        nap += 1
        if nap > napmax:
            logic = 4
            fn = fg
            xn = xn + tg * d
            return x, fn, g, nap, logic
        indic = 4
        #
        #     --- appel simulateur
        #
        f, g = simulator.simul(x, run_id=nap, **kwargs)
        # , 1999 + nap - 1, io

        indic = 1

        if indic == 0:
            #
            #         --- arret demande par l'utilisateur
            #
            logic = 5
            fn = f
            xn = x
            return x, fn, g, nap, logic

        if indic < 0:
            #
            #         --- les calculs n'ont pas pu etre effectues par le
            # simulateur
            #
            td = t
            indicd = indic
            logic = 0
            if imp >= 4:
                info("     mlis0  t=" + str(t) + " indic=" + str(indic), io)
            t = tg + 0.1 * (td - tg)

        else:
            #
            #     --- les tests elementaires sont faits, on y va
            #
            ps = np.dot(d, g)
            fp = ps
            #
            #     --- premier test de Wolfe
            #
            ffn = f - fn
            goto900 = 0
            if ffn > t * tesf:
                td = t
                fd = f
                fpd = fp
                indicd = indic
                logic = 0
                if imp >= 4:
                    info(
                        "     mlis0  t="
                        + str(t)
                        + " ffn="
                        + str(ffn)
                        + " fp="
                        + str(fp),
                        io,
                    )

            else:
                #
                #     --- test 1 ok, donc deuxieme test de Wolfe
                #
                if imp >= 4:
                    info(
                        "     mlis0  t="
                        + str(t)
                        + " ffn="
                        + str(ffn)
                        + " fp="
                        + str(fp),
                        io,
                    )

                if fp > tesd:
                    logic = 0
                    fn = f
                    xn = x
                    return x, fn, g, nap, logic

                if logic != 0:
                    #
                    #     --- test 2 ok, donc pas serieux, on sort
                    #
                    fn = f
                    xn = x
                    return x, fn, g, nap, logic

                tg = t
                fg = f
                fpg = fp
                if not td:
                    #
                    #              extrapolation
                    #
                    taa = t
                    gauche = (1.0 + barmin) * t
                    droite = 10.0 * t
                    t = ecube(t, f, fp, ta, fa, fpa, gauche, droite)
                    ta = taa
                    goto900 = 1
                    if t >= tmax:
                        logic = 1
                        t = tmax
            #
            #              interpolation
            #
            if not goto900 and indica <= 0:
                ta = t
                t = 0.9 * tg + 0.1 * td
                goto900 = 1

            if not goto900:
                test = barr * (td - tg)
                gauche = tg + test
                droite = td - test
                taa = t
                t = ecube(t, f, fp, ta, fa, fpa, gauche, droite)
                ta = taa
                if gauche < t < droite:
                    barr = max(barmin, barr / barmul)
                else:
                    barr = min(barmul * barr, barmax)
            #
            # --- fin de boucle
            #     - t peut etre bloque sur tmax
            #       (venant de l'extrapolation avec logic=1)
            #
            fa = f
            fpa = fp

        indica = indic
        #
        # --- faut-il continuer ?
        #
        if td != 0.0:
            goto950 = 0
            if (td - tg) >= tmin:
                #
                #     --- limite de precision machine (arret de secours) ?
                #
                for i in range(n):
                    z = xn[i] + t * d[i]
                    if z != xn[i] and z != x[i]:
                        goto950 = 1
                        break
            if not goto950:
                #
                # --- arret sur dxmin ou de secours
                #
                logic = 6
                #
                #     si indicd<0, derniers calculs non faits par simul
                #
                if indicd < 0:
                    logic = indicd
                #
                #     si tg=0, xn = xn_depart,
                #     sinon on prend xn=x_gauche qui fait decroitre f
                #
                if tg != 0.0:
                    fn = fg
                    xn = xn + tg * d
                if imp <= 0:
                    return x, fn, g, nap, logic

                info(
                    "  mlis0   stop_on_tmin   step   functions   derivatives",
                    io,
                )
                info(
                    "  mlis0   " + str(tg) + "  " + str(fg) + " " + str(fpg),
                    io,
                )

                if logic == 6:
                    info(
                        "  mlis0   "
                        + str(td)
                        + "  "
                        + str(fg)
                        + " "
                        + str(fpd),
                        io,
                    )

                if logic == 7:
                    info(
                        "  mlis0   " + str(td) + "    indic=" + str(indicd), io
                    )

                return x, fn, g, nap, logic
        #
        #               recopiage de x et boucle
        #
        x = xn + t * d

    return x, fn, g, nap, logic


#
# -----------------------------------------------------------------------
#


def descentdir(
    n, sscale, nm, depl, jmin, jmax, precos, diag, alpha, ybar, sbar
):
    """
        Decrease descent
       calcule le produit H.g ou
           . H est une matrice construite par la formule de bfgs inverse
             a nm memoires a partir de la matrice diagonale diag
             dans un espace hilbertien dont le produit scalaire
             est donne par prosca
             (cf. J. Nocedal, Math. of Comp. 35/151 (1980) 773-782)
           . g est un vecteur de dimension n (en general le gradient)

       la matrice diag apparait donc comme un preconditionneur diagonal

       depl = g (en entree), = H g (en sortie)

       la matrice H est memorisee par les vecteurs des tableaux
       ybar, sbar et les pointeurs jmin, jmax

       alpha(nm) est une zone de travail

       izs(1),rzs(1),dzs(1) sont des zones de travail pour prosca

           arguments

       logical sscale
       integer n,nm,jmin,jmax
       real depl(n),precos,diag(n),alpha(nm),ybar(n,10),sbar(n,10),
      &     aux(n)
       external prosca,ctonb,ctcab

           variables locales

       integer jfin,i,j,jp
       real r
       double precision ps

    """
    jfin = jmax
    if jfin < jmin:
        jfin = jmax + nm
    #
    #         phase de descente
    #
    for j in range(jfin, jmin - 1, -1):
        jp = j
        if jp >= nm:
            jp -= nm
        ps = np.dot(depl, sbar[jp, :])
        alpha[jp] = ps
        depl = depl - ps * ybar[jp, :]
    #
    #         preconditionnement
    #
    if sscale:
        depl = depl * precos
    else:
        depl = depl * diag
    #
    #         remontee
    #
    for j in range(jmin, jfin + 1):
        jp = j
        if jp >= nm:
            jp = jp - nm

        ps = np.dot(depl, ybar[jp, :])
        r = alpha[jp] - ps
        depl = depl + r * sbar[jp, :]

    return depl


def ecube(t, f, fp, ta, fa, fpa, tlower, tupper):
    """
    --- arguments

       real t,f,fp,ta,fa,fpa,tlower,tupper

    --- variables locales

       real sign,den,anum
       double precision z1,b,discri

             Using f and fp at t and ta, computes new t by cubic formula
             safeguarded inside [tlower,tupper].
    """

    z1 = fp + fpa - 3.0 * (fa - f) / (ta - t)
    b = z1 + fp
    #
    #              first compute the discriminant (without overflow)
    #
    if abs(z1) <= 1.0:
        discri = z1 * z1 - fp * fpa
    else:
        discri = fp / z1
        discri = discri * fpa
        discri = z1 - discri
        if z1 * discri > 0.0:
            discri = z1 * discri
        else:
            discri = -1.0

    if discri < 0.0:
        if fp < 0.0:
            t = tupper
        if fp >= 0.0:
            t = tlower

    else:
        #
        #  discriminant nonnegative, compute solution (without overflow)
        #
        discri = np.sqrt(discri)
        if t - ta < 0.0:
            discri = -discri

        sign = (t - ta) / abs(t - ta)
        if b * sign > 0.0:
            t = t + fp * (ta - t) / (b + discri)

        else:
            den = z1 + b + fpa
            anum = b - discri
            if abs((t - ta) * anum) < (tupper - tlower) * abs(den):
                t = t + anum * (ta - t) / den
            else:
                t = tupper

    t = max(t, tlower)
    t = min(t, tupper)

    return t
