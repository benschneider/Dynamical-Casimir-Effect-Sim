from scipy.optimize import curve_fit, minimize
from scipy import constants as cons
import numpy as np
from numpy import pi, cos, sin
from PyGnuplot import gp
from struct import pack, unpack
# import numexpr as ne

class emptyclass1():
    pass


fit_variables = emptyclass1()

def load_data(fileName):
    a = open(fileName, 'r+')
    b = a.readlines()
    data = np.zeros([2, len(b[:-3])])
    for i, e in enumerate(b[:-3]):
        strings = e.split(' ')
        data[0, i] = (float(strings[0]))  # raw xaxis
        data[1, i] = (float(strings[1]))  # raw yaxis
    return data

def get_s11(flux, Ic, Cap, R=1e99):
    flux0 = cons.h/(2*cons.e)
    f0 = 4.1e9
    z0 = 50
    L = flux0 / (Ic * 2.0 * pi * np.abs(np.cos(pi * flux)) + 1e-90)
    Ysq = (1.0 / R + 1.0 / (1j * 2.0 * pi * f0 * L + 1j * 1e-90) + 1j * 2.0 * pi * f0 * Cap)
    zsq = 1.0 / Ysq
    s11 = (zsq - z0) / (zsq + z0)
    return s11  # returns the complex s11 response

def fitFunc_ang(flux, Ic, Cap, offset1=0.0, slope=0.0):
    s11 = get_s11(flux, Ic, Cap)
    return np.angle(s11) + offset1 - slope * flux

def get_fitvals(obj1=fit_variables, file_ang = '1169S11_p_4p1.dat'):
    data = load_data(file_ang)
    R = 3000.0
    Ic = 2.4e-6
    Cap = 6e-13
    flux = data[0]
    offset = 0.0
    slope = 0.1
    scale2 = 1.0  # for the fitting
    flux0 = cons.h/(cons.e*2)
    iguess = [Ic, Cap, offset, slope]
    popt, pcov = curve_fit(fitFunc_ang, flux, data[1], p0=iguess)  # ,
    obj1.Ic = popt[0]
    obj1.Cap = popt[1]
    obj1.off = popt[2]
    obj1.slope = popt[3]
    return popt, data[1], data[0]

def get_ang_resp(x, obj=fit_variables):
    # uses fitted values to return an anglular value
    Ic = obj.Ic
    Cap = obj.Cap
    off = obj.off
    slope = obj.slope
    return fitFunc_ang(x, Ic, Cap, off, slope)

if __name__ == "__main__":
    [a,b,c,d], ydata, flux = get_fitvals(file_ang = '1169S11_p_4p1.dat')
    yfit = get_ang_resp(flux) # get angular response
    fig1 = gp()
    fig1.s([flux, yfit, ydata], "fitplot.dat")
    fig1.c('plot "fitplot.dat" u 1:2 w l t "fit"')
    fig1.c('replot "fitplot.dat" u 1:3 w l t "phase"')
