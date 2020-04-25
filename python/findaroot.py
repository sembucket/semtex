#!/usr/bin/env python
#
# $Id: findaroot.py,v 9.1 2019/05/30 06:36:10 hmb Exp $
# ----------------------------------------------------------------------------
import math


# ============================================================================
def secant (f, x0, xstep=0.0, ytol=1.0e-6, maxit=30):
    '''
    # Secant method for finding root of f(x).
    # Adapted from stackoverflow answer by Hugh Bothwell.
    #
    # f:     a function that takes a single argument.
    # x0:    an initial guess for the root.
    # xstep: if xstep is zero our initial step will be 2% of x0,
    #        otherwise we use the suggested xstep.
    # tol:   the convergence tolerance on f(x).
    # maxit: the maximum allowed number of iterations.
    '''

    if (xstep == 0.0): x1 = 1.02 * x0
    else:              x1 = x0 + xstep

    fx0 = f(x0)

    if (abs (fx0) < ytol): return x0
    
    fx1 = f(x1)

    for _ in range (maxit):
        if (abs (fx1) < ytol): return x1
        x2 = (x0 * fx1 - x1 * fx0) / (fx1 - fx0)
        x0,  x1  = x1,  x2
        fx0, fx1 = fx1, f(x2)
    raise Exception ("secant exceeded maximum iterations")


# ============================================================================
def zriddr (f, x1, x2, xtol=1.0e-6, ytol=1.0e-6, maxit=30):
    '''
    # Ridder's method for finding root of f(x) adapted from Numerical Recipes.
    # x1 and x2 must bracket root.
    '''
    fl = f(x1)
    fh = f(x2)
    if (abs(fl) <= ytol): return x1
    if (abs(fh) <= ytol): return x2    
    if ( ((fl > 0.0) and (fh < 0.0)) or ((fl < 0.0) and (fh > 0.0)) ):
        xl = x1
        xh = x2
        ans = -9.99e99
        for _ in range (maxit):
            xm = 0.5*(xl+xh)
            fm = f(xm)
            if (abs(fm) <= ytol): return xm            
            # print "xm = ", xm, "f(xm) = ", fm            
            s = math.sqrt(fm*fm-fl*fh)
            if (s == 0.0): return ans
            xnew = xm+(xm-xl)*((1.0 if (fl >= fh) else -1.0) * fm/s)
            if (abs(xnew-ans) <= xtol): return ans
            ans = xnew
            fnew = f(ans)
            if (fnew == 0.0): return ans
            if (math.copysign(fm,fnew) != fm):
                xl = xm
                fl = fm
                xh = ans
                fh = fnew
            elif (math.copysign(fl,fnew) != fl):
                xh = ans
                fh = fnew
            elif (math.copysign(fh,fnew) != fh):
                xl = ans
                fl = fnew
            else:
                raise Exception ("never get here")
            if (abs(xh-xl) <= xtol): return ans
        raise Exception ("zriddr exceed maximum iterations")
    else:
        if (fl == 0.0): return xl
        if (fh == 0.0): return xh
        raise Exception ("root must be bracketted in zriddr")
   

# ============================================================================
def zbrak (f, x1, x2, ytol=1e-6, maxit=3):
    '''
    # Geometric expansion of initial guessed range until root
    # is bracketted by returned values x1 and x2.
    # F(x1) and f(x2) are also returned.
    # Adapted from Numerical Recipes.
    '''
    if (x1 == x2) : raise Exception ("bad initial range in zbrak")
    f1 = f(x1)
    f2 = f(x2)
    if ((abs(f1) <= ytol) or (abs(f2) <= ytol)) : return x1, x2, f1, f2
    for i in range (maxit):
        if (f1*f2 <= 0.0): return x1, x2, f1, f2
        if (abs(f1) < abs(f2)):
            tmp = x1
            x1 += 1.6*(x1-x2)
            x2 = tmp
            f2 = f1
            f1 = f(x1)
        else:
            tmp = x2
            x2 += 1.6*(x2-x1)
            x1 = tmp
            f1 = f2
            f2 = f(x2)
    raise Exception ("too many iterations in zbrak")


# ============================================================================
def zriddrbrak (f, x1guess, x2guess, xtol=1.0e-6, ytol=1.0e-6, maxit=30):
    '''
    # Ridder's method for finding root of f(x) adapted from Numerical Recipes.
    # x1guess and x2guess are used as initial input to bracketting code.
    '''
    x1, x2, fl, fh = zbrak (f, x1guess, x2guess, ytol)
    if (abs(fl)    <= ytol): return x1
    if (abs(fh)    <= ytol): return x2
    if (abs(x1-x2) <= xtol): return 0.5*(x1+x2)
    if ( ((fl > 0.0) and (fh < 0.0)) or ((fl < 0.0) and (fh > 0.0)) ):
        xl = x1
        xh = x2
        ans = -9.99e99
        for i in range (maxit):
            xm = 0.5*(xl+xh)
            fm = f(xm)
            if (abs(fm) <= ytol): return xm
#            print "xm = ", xm, "f(xm) = ", fm            
            s = math.sqrt(fm*fm-fl*fh)
            if (s == 0.0): return ans
            xnew = xm+(xm-xl)*((1.0 if (fl >= fh) else -1.0) * fm/s)
            if (abs(xnew-ans) <= xtol): return ans
            ans = xnew
            fnew = f(ans)
            if (abs(fnew) <= ytol): return ans
            if (math.copysign(fm,fnew) != fm):
                xl = xm
                fl = fm
                xh = ans
                fh = fnew
            elif (math.copysign(fl,fnew) != fl):
                xh = ans
                fh = fnew
            elif (math.copysign(fh,fnew) != fh):
                xl = ans
                fl = fnew
            else:
                raise Exception ("never get here")
            if (abs(xh-xl) <= xtol): return ans
        raise Exception ("zriddrbrak exceed maximum iterations")
    else:
        if (fl == 0.0): return xl
        if (fh == 0.0): return xh
        raise Exception ("root must be bracketted in zriddrbrak")


# ============================================================================
# Testing, testing...
def testfun (x):
    # This has one root at x = 0.56714329.
    return x - math.exp(-x)

#root = secant (testfun, 20)
#print "root using secant method: ", root
#
root = zriddr (testfun, 0, 20)
print "root using Ridder's method: ", root
#
#x1, x2, f1, f2 = zbrak (testfun, 0.2, 0.3)
#print "bracketting interval returned by zbrak: ", x1, x2
#
root = zriddrbrak (testfun, 0.2, 0.4, 1e-4, 1e-4)
print "root using Ridder's method with bracketting: ", root

