from scipy.optimize import curve_fit, minimize
from scipy import constants as cons
import numpy as np
from numpy import pi, cos, sin
import PyGnuplot as gp

import fit_mag_phase as fit_data
import make_hoa as hoa
from loadsavemtx import savemtx, loadmtx
'''
a
fit_mag_phase.py
fit phase response, have continous dynamic phase data ready from fit

b
make_hoa.py:
    Multiply signal with phase data -> do FFT on resulting response
 save higher order amplitudes to a file (hoa) for later usage.
(hoa file: contains hoa values for flux and amplitude values)

c
This script:
multiply hoa (higher order amplitudes) with equation to obtain responses: make_hoa.py
 have 4 different axes: detector frequency, flux DC, flux AC, plasmafreq
 Results of DCE photon rates are then separated into 2 maps,
    1. v.s. plasmafrequency for a fixed detector frequency (E.g. 4.1GHz),
    2. v.s. detector frequency for a fixed plasma (E.g. 14.5Ghz)

'''


def LJDC(fluxdc, Ic=3.4e-6):
    # fluxac and fluxdc are normalized by flux0
    flux0 = cons.h/(2*cons.e)
    Ej_dc = Ic*flux0/(2*np.pi)
    return (flux0/(2*pi))**2 *1/(Ej_dc*np.abs(cos(pi*fluxdc))) # normalized fluxdc

def get_beta(D, w, wd=8.9e9, wp=40e9, fluxdc=0.475, Ic=3.4e-6, F0=np.pi, dfdc=0):  #, Ic=3.4e-6):
    # F0 is the conversion factor from hoa phase to an effective F0*D=fluxac
    flux0 = cons.h/(2*cons.e)
    Ldc = LJDC(fluxdc, Ic)  # SQUID inductance
    Z0 = 50
    part1 = F0*D*2j*pi*D # fluxpump strength
    part2 = np.sqrt(w*(wd-w))/(wd*Z0)
    part3 = 1-(w**2)/(wp**2)+1j*Ldc*w/Z0  # SQUID reflection eqn r(w)(negative w)
    part4 = 1-(wd-w)**2/(wp**2)-1j*Ldc*(wd-w)/Z0  # SQUID reflection eqn r(w)(positive w)
    beta = part1*part2/(part3*part4)
    return np.nan_to_num(beta)

def get_response_vs_pfreq(w, wd, pfreq, dc_flux, ac_flux, harmonics):
    phot_map = 1j*np.zeros([len(dc_flux),len(ac_flux),len(pfreq)])
    hoa_str = np.zeros(len(ac_flux))
    for i, fluxdc in enumerate(dc_flux):
        print(i)
        for k, wp in enumerate(pfreq):
            for harmonic in range(harmonics):
                wd_eff = wd + wd*harmonic
                hoa_str = hoa_matrix[i, :, harmonic]  # get vector
                beta = get_beta(hoa_str, w, wd_eff, wp, fluxdc)
                phot_map[i,:,k] += beta
    return np.abs(phot_map)**2


def get_response_vs_w(ws, wd, pfreq, dc_flux, ac_flux, harmonics):
    phot_map = 1j*np.zeros([len(dc_flux),len(ac_flux),len(ws)])
    hoa_str = np.zeros(len(ac_flux))
    for i, fluxdc in enumerate(dc_flux):
        print(i)
        for k, w in enumerate(ws):
            wp = pfreq[0]
            for harmonic in range(harmonics):
                wd_eff = wd + wd*harmonic
                hoa_str = hoa_matrix[i, :, harmonic]  # get vector
                beta = get_beta(hoa_str, w, wd_eff, wp, fluxdc)
                phot_map[i,:,k] += beta
    return np.abs(phot_map)**2

fit_data.get_fitvals()

# Uncomment this code for creation of 2D maps
dc_flux = np.linspace(-0.75, 0.70, 291)
ac_flux = np.linspace(0, 1.0, 1001)
harmonics = 20  # only use first 20 harmonics
wd = 8.9e9
timeaxis = np.linspace(0, 6e-9, 20001)  # for FFT
w = 4.8e9  # detector frequency
ws = np.linspace(0.01e9, 15.01e9, 151)
pfreq = [14.5e9] # for single plasmafrequency map
pfreqs = np.linspace(1, 31, 121)*1e9  # plasmafrequencies
output_file1 = 'output/mtx/hoa_matrix_wide_20.mtx'
output_file2 = 'output/mtx/dce_p_20.mtx'
output_file3 = 'output/mtx/dce_w_20.mtx'

# uncomment this to create a higher harmonics map
hoa_matrix = hoa.make_response(ac_flux, dc_flux, wd, timeaxis, fit_data, harmonics=harmonics, plot=False)
h = ('Units,ufo,Harmonics,' + str(0) + ',' + str(harmonics) +
     ',Pump,' + str(ac_flux[0]) + ',' + str(ac_flux[-1]) +
     ',FluxPos,' + str(dc_flux[0]) + ','+str(dc_flux[-1])+'')
savemtx(output_file1, hoa_matrix, h)

# loading and using a generated higher harmonics map vs plasma frequency
hoa_matrix, h = loadmtx(output_file1)
print('calculating beta map with '+str(harmonics)+' harmonics.')
print('and '+str(len(pfreq))+' Plasma-frequencies.')
DCE_map = get_response_vs_pfreq(w, wd, pfreqs, dc_flux, ac_flux, harmonics)
header2 = ('Units,Photon,Plasma Freq,' + str(pfreqs[0]) + ',' + str(pfreqs[-1]) +
          ',Pump,' + str(ac_flux[0]) + ',' + str(ac_flux[-1]) +
          ',FluxPos,' + str(dc_flux[0]) + ','+str(dc_flux[-1])+'')
savemtx(output_file2, DCE_map, header2)


# loading and using a generated higher harmonics map vs detector frequency
print('calculating ws map')
hoa_matrix, h = loadmtx(output_file1)
DCE_map_ws = get_response_vs_w(ws, wd, pfreq, dc_flux, ac_flux, harmonics)
header3 = ('Units,Photon,Det. Freq,' + str(ws[0]) + ',' + str(ws[-1]) +
          ',Pump,' + str(ac_flux[0]) + ',' + str(ac_flux[-1]) +
          ',FluxPos,' + str(dc_flux[0]) + ','+str(dc_flux[-1])+'')
savemtx(output_file3, DCE_map_ws, header3)

