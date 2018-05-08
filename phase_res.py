#! /Library/Frameworks/Python.framework/Versions/3.6/bin/python3 -i
from scipy.optimize import curve_fit, minimize
import numpy as np
from numpy import pi, cos, sin
# import matplotlib.pyplot as plt
import PyGnuplot as gp
from struct import pack, unpack

plt.ion()


def savemtx(filename, data, header='Units,ufo,d1,0,1,d2,0,1,d3,0,1'):
    with open(filename, 'wb') as f:
        f.write(bytes(header + '\n', 'utf-8'))
        mtxshape = data.shape
        line = str(mtxshape[2]) + ' ' + str(mtxshape[1]) + ' ' + str(mtxshape[0]) + ' ' + '8'
        f.write(bytes(line + '\n', 'utf-8'))  # 'x y z 8 \n'
        for ii in range(mtxshape[2]):
            for jj in range(mtxshape[1]):
                content = pack('%sd' % mtxshape[0], *data[:, jj, ii])
                f.write(content)
        f.close()


def load_data(fileName):
    a = open(fileName, 'r+')
    b = a.readlines()
    data = np.zeros([3, len(b[:-3])])
    for i, e in enumerate(b[:-3]):
        strings = e.split(' ')
        data[0, i] = (float(strings[0]))  # raw xaxis
        data[2, i] = (float(strings[1]))  # raw yaxis
    return data


def return_2p_value(xaxis, ydata, xvalin):
    # approximate intermediate value
    # if (xval < xaxis[1]) or (xval > xaxis[-2]):
    xval = (xvalin - 0.3) % 1.0 - 0.7  # wrap values into: (-0.3 to 0.7)
    xlength = len(xaxis)
    xstep = (xaxis[-1] - xaxis[0]) / (xlength-1)
    idx = (xval-xaxis[0])/xstep
    int_idx = int(idx)
    resid_idx = idx - int_idx
    resid = (ydata[int_idx+1] - ydata[int_idx]) * resid_idx
    return ydata[int_idx] + resid


def get_s11(flux, Ic, Cap, R=1e99):
    flux0 = 2.07e-15
    f0 = 4.1e9
    z0 = 50
    L = flux0 / (Ic * 2.0 * pi * np.abs(np.cos(pi * flux)) + 1e-90)
    Ysq = (1.0 / R + 1.0 / (1j * 2.0 * pi * f0 * L + 1j * 1e-90) + 1j * 2.0 * pi * f0 * Cap)
    zsq = 1.0 / Ysq
    s11 = (zsq - z0) / (zsq + z0)
    return s11


def fitFunc_mag(flux, Ic, Cap, R, offset2, scales2):
    s11 = get_s11(flux, Ic, Cap, R)
    return scales2 * np.abs(s11) + offset2


def fitFunc_ang(flux, Ic, Cap, offset1=0.0, slope=0.0):
    s11 = get_s11(flux, Ic, Cap)
    return np.angle(s11) + offset1 - slope * flux


def parabola(a, f, x):
    if a < 1e-5:
        return x * 0.0
    elif f < 1e-20:
        return x * 0.0
    else:
        b = f / 2.0
        z = (-a / b**2 * (x - b)**2 + a)
        z[z < 0] = 0
    return z


def get_full_parabola(freqaxis, fft_drive, scale=0.01):
    parabola_full = np.zeros_like(freqaxis)
    for i, amp in enumerate(fft_drive):
        freq = freqaxis[i]
        parabola_full += parabola(amp * scale, freq, freqaxis)
    return parabola_full


def get_drive(amplitude, dc_offset, omega0, timeaxis):
    drive = np.zeros([4, resolution])
    for i, jj in enumerate(timeaxis):
        signal = amplitude * np.sin(jj * omega0) + dc_offset
        drive[0, i] = signal
        # drive[1, i] = return_2p_value(flux, data_mag[2], signal)
        # drive[3, i] = return_2p_value(flux, data[2], signal)
        drive[1, i] = fitFunc_mag(signal, popt2[0], popt2[1], popt2[2], popt2[3], popt2[4])
        drive[3, i] = fitFunc_ang(signal, Ic, Cap, offset, slope)
    return drive


def get_fft_responses(drive):
    loss = np.mean(drive[1])
    m0 = np.min(drive[0]) + (np.max(drive[0]) - np.min(drive[0])) / 2.0
    m3 = np.min(drive[3]) + (np.max(drive[3]) - np.min(drive[3])) / 2.0
    fft_signal = np.abs(np.fft.rfft(drive[0] - m0, norm="ortho"))
    fft_mirror = np.abs(np.fft.rfft(drive[3] - m3, norm="ortho")) * loss
    # fft_loss = np.abs(np.fft.rfft(drive[1], norm="ortho"))
    return fft_signal, fft_mirror


def get_LogNegNum(CovM_in):
    CovM = CovM_in * 2
    V = np.linalg.det(CovM)
    A = np.linalg.det(CovM[:2, :2])
    B = np.linalg.det(CovM[2:, 2:])
    C = np.linalg.det(CovM[:2, 2:])
    sigma = A + B - 2.0 * C
    vn = np.sqrt(sigma / 2.0 - np.sqrt(sigma * sigma - 4.0 * V) / 2.0)
    if C == 0:
        return 0.0
    else:
        return -np.log(2.0 * vn) if (np.log(2.0 * vn)) < 0.0 else 0.0


def get_Negativity(n12, nth=0.0):
    M_test = np.array([[(0.25+(n12+nth)/2.0), 0, np.sqrt(n12)/2.0, 0],
                       [0, (0.25 + (n12+nth)/2.0), 0, -np.sqrt(n12)/2.0],
                       [np.sqrt(n12)/2.0, 0, (0.25+(n12+nth)/2.0), 0],
                       [0, -np.sqrt(n12)/2.0, 0, (0.25+(n12+nth) / 2.0)]])
    N = get_LogNegNum(M_test)
    return N

fileName_ang = '1157_S11_Avg_phase.linecut'
fileName_mag = '1157_S11_Avg_mag.linecut'
output_file1 = 'sim_all_parabolas_2.mtx'
output_file2 = 'sim_pure_parabola_2.mtx'
output_file3 = 'negativity_fit_2.mtx'
output_file4 = 'TMS_NONTMS_ratio.mtx'

data = load_data(fileName_ang)
data_mag = load_data(fileName_mag)
x0 = -0.792
x1 = 4.03
data[1] = ((data[0] - x0) / (x1 - x0)) - 0.5  # corrected xaxis
# Fitting Squid Phase response
R = 300.0
Ic = 2.4e-6
Cap = 3e-13
flux0 = 2.07e-15    # Tm^2; Flux quanta: flux0 =  h / (2*charging energy)
flux = np.linspace(-0.75, 0.7011, 701)  # * flux0
offset = -4.5
slope = 0.1
offset2 = 0.0
scale2 = 1.0
iguess = [Ic, Cap, offset, slope]
iguess2 = [Ic, Cap, R, offset2, scale2]
popt, pcov = curve_fit(fitFunc_ang, flux, data[2], p0=iguess)
data_mag[2] = data_mag[2] / np.max(data_mag[2])
popt2, pcov2 = curve_fit(fitFunc_mag, flux, data_mag[2], p0=iguess2, maxfev=500000, xtol=1e-36)
# Only use the fit to obtain intermediate and lower noise data
# print(popt)
# print(pcov)
Ic, Cap, offset, slope = popt
print(popt2)
print(pcov2)
resolution = 2**11
pumpfreq = 8.9e9
omega0 = 2.0 * np.pi * pumpfreq
timeaxis = np.linspace(0.0, 5.0e-9, resolution)
freqaxis = np.fft.rfftfreq(timeaxis.shape[-1], (timeaxis[1] - timeaxis[0]))
#  gap_freq = 88e9  # kept the FFT nyquist limit below this
pump_idx = np.argmin(abs(freqaxis - pumpfreq))
scale = 0.01
fluxpoints = 561
powerpoints = 101
dc_offsets = np.linspace(-0.65, 0.75, fluxpoints)
amplitudes = np.linspace(0.0, 0.1, powerpoints)
parabolas = np.zeros([fluxpoints, powerpoints, len(freqaxis)])
parabolas2 = np.zeros([fluxpoints, powerpoints, len(freqaxis)])
negativity = np.zeros([fluxpoints, powerpoints, len(freqaxis)])
# TMS_ratio = np.zeros([fluxpoints, powerpoints, len(freqaxis)])

for kk, dc_offset in enumerate(dc_offsets):
    print(kk)
    for jj, amplitude in enumerate(amplitudes):
        drive = get_drive(amplitude, dc_offset, omega0, timeaxis)
        fft_signal, fft_mirror = get_fft_responses(drive)
        amp_d = fft_mirror[pump_idx]  # FFT amp component at drive
        # Flux parabola with all Frequencies
        parabolas[kk, jj] = get_full_parabola(freqaxis, fft_mirror, scale)
        # Flux parabola only at pump
        parabolas2[kk, jj] = parabola(amp_d * scale, pumpfreq, freqaxis)
        # parabolas[2, jj] = (parabolas[1, jj] + 1e-199) / (parabolas[0, jj] + 1e-90)
        # parabolas[3, jj] = parabolas[1, jj] - parabolas[0, jj]
        added_photons = parabolas[kk, jj] - parabolas2[kk, jj]
        for ii, ph in enumerate(parabolas2[kk, jj]):
            negativity[kk, jj, ii] = get_Negativity(ph, nth=added_photons[ii])
            # TMS_ratio[kk, jj, ii] = parabolas2[kk, jj]/parabolas[kk, jj]

header = ('Units,ufo,Frequency,' + str(freqaxis[0]) + ',' + str(freqaxis[-1]) +
          ',Pump,' + str(amplitudes[0]) + ',' + str(amplitudes[-1]) +
          ',FluxPos,' + str(dc_offsets[0]) + ','+str(dc_offsets[-1])+'')
# savemtx(output_file1, parabolas, header)
# savemtx(output_file2, parabolas2, header)

added_photons = parabolas - parabolas2
# shape = parabolas2.shape
# for kk in range(parabolas2.shape[0]):
#     print(kk)
#     for jj in range(parabolas2.shape[1]):
#         for ii, ph in enumerate(parabolas2[kk, jj]):
#             negativity[kk, jj, ii] = get_Negativity(ph, nth=added_photons[kk, jj, ii])

# savemtx(output_file3, negativity, header)
TMS_ratio = parabolas2 / parabolas
savemtx(output_file4, TMS_ratio, header)

m0 = np.min(drive[0]) + (np.max(drive[0]) - np.min(drive[0])) / 2.0
m3 = np.min(drive[3]) + (np.max(drive[3]) - np.min(drive[3])) / 2.0
plt.figure(1)
plt.plot(flux, data_mag[2])
plt.plot(flux, fitFunc_mag(flux, popt2[0], popt2[1], popt2[2], popt2[3], popt2[4]))
plt.figure(2)
plt.plot(flux, data[2])
plt.plot(flux, fitFunc_ang(flux, popt[0], popt[1], popt[2], popt[3]))
plt.figure(4)
plt.plot(timeaxis, drive[0] - m0)
plt.figure(5)
plt.plot(timeaxis, drive[3] - m3)
plt.figure(6)
plt.plot(freqaxis, fft_signal)
plt.figure(7)
plt.plot(freqaxis, fft_mirror)
plt.figure(8)
plt.plot(freqaxis, parabolas[0, jj])
