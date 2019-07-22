from scipy.optimize import curve_fit, minimize
from scipy import constants as cons
import numpy as np
from numpy import pi, cos, sin
import PyGnuplot as gp
from struct import pack, unpack
# import numexpr as ne

class emptyclass1():
    pass


fit_variables = emptyclass1()

def load_data(fileName):
    a = open(fileName, 'r+')
    b = a.readlines()
    data = np.zeros([3, len(b[:-3])])
    for i, e in enumerate(b[:-3]):
        strings = e.split(' ')
        data[0, i] = (float(strings[0]))  # raw xaxis
        data[2, i] = (float(strings[1]))  # raw yaxis
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

# def get_s11_numexpr(flux, Ic, Cap, R=1e99):
#     # 4ms faster per execution
#     flux0 = cons.h/(2*cons.e)
#     f0 = 4.1e9
#     z0 = 50
#     zsq = ne.evaluate("1.0 / ((1.0 / R + 1.0 / (1j * 2.0 * pi * f0 *(flux0 / (Ic * 2.0 * pi * abs(cos(pi * flux)) + 1e-90)) + 1j * 1e-90) + 1j * 2.0 * pi * f0 * Cap))")
#     s11 = (zsq - z0) / (zsq + z0)
#     return s11  # returns the complex s11 response

def fitFunc_mag(flux, Ic, Cap, R, offset2, scales2):
    s11 = get_s11(flux, Ic, Cap, R)
    return scales2 * np.abs(s11) + offset2

def fitFunc_ang(flux, Ic, Cap, offset1=0.0, slope=0.0):
    s11 = get_s11(flux, Ic, Cap)
    return np.angle(s11) + offset1 - slope * flux

def get_fitvals(obj1=fit_variables):
    fileName_ang = 'SQUID_phase.dat'
    fileName_mag = 'SQUID_magnitude.dat'
    data = load_data(fileName_ang)
    data_mag = load_data(fileName_mag)
    x0 = -0.792
    x1 = 4.03
    data[1] = ((data[0] - x0) / (x1 - x0)) - 0.5  # corrected xaxis
    # Fitting Squid Phase response
    R = 3000.0
    Ic = 2.4e-6
    Cap = 6e-13
    flux0 = cons.h/(cons.e*2)
    flux = np.linspace(-0.75, 0.7011, 701)  # * flux0
    offset = -4.5
    slope = 0.1
    offset2 = 0.0
    scale2 = 1.0  # for the fitting
    iguess = [Ic, Cap, offset, slope]
    iguess2 = [Ic, Cap, R, offset2, scale2]
    popt, pcov = curve_fit(fitFunc_ang, flux, data[2], p0=iguess)  # ,
                           #bounds=([1e-9,1e-16,-10,-0.5],[1e-4, 1e-9,10,0.5]),
                           #maxfev=5000000, xtol=1e-10)
    # data_mag[2] = data_mag[2] / np.max(data_mag[2])
    # popt2, pcov2 = curve_fit(fitFunc_mag, flux, data_mag[2], p0=iguess2, maxfev=500000, xtol=1e-36)
    # Ic, Cap, offset, slope = popt
    obj1.Ic = popt[0]
    obj1.Cap = popt[1]
    obj1.off = popt[2]
    obj1.slope = popt[3]
    return popt, data[2]
    #print(popt2)
    #print(pcov2)

def get_ang_resp(x, obj=fit_variables):
    # uses fitted values to return an anglular value
    Ic = obj.Ic
    Cap = obj.Cap
    off = obj.off
    slope = obj.slope
    return fitFunc_ang(x, Ic, Cap, off, slope)

if __name__ == "__main__":
    [a,b,c,d], curvepoints = get_fitvals()
    x = np.linspace(-0.75,0.7011,701)
    y = get_ang_resp(x) # get angular response
    gp.s([x,y, curvepoints], "fitplot.dat")
    gp.c('plot "fitplot.dat" u 1:2 w l t "fit"')
    gp.c('replot "fitplot.dat" u 1:3 w l t "phase"')
