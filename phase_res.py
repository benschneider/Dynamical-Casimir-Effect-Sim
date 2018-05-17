from scipy.optimize import curve_fit, minimize
import numpy as np
from numpy import pi, cos, sin
import PyGnuplot as gp
from struct import pack, unpack

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

def get_first_parabola(freqaxis, fft_drive, scale=0.01, fcent=8.9e9, cutoff=0):
    # takes only values around the first pump
    parabola_full = np.zeros_like(freqaxis)
    for i, amp in enumerate(fft_drive):
        freq = freqaxis[i]
        if (amp > cutoff) and (freq < fcent+1e9) and (freq > fcent-1e9):
            # taking a witdh around fcent to minimize fft artefacts
            parabola_full += parabola(amp * scale, freq, freqaxis)
    return parabola_full

def get_full_parabola(freqaxis, fft_drive, scale=0.01, cutoff=0):
    parabola_full = np.zeros_like(freqaxis)
    for i, amp in enumerate(fft_drive):
        freq = freqaxis[i]
        if (amp > cutoff) and (freq > 1e9):
            parabola_full += parabola(amp * scale, freq, freqaxis)
    return parabola_full


def get_drive(amplitude, dc_offset, omega0, timeaxis):
    # 0: Input amplitude v.s. time
    # 1: Corresponding magnitude v.s. time
    # 2: Corresponding phase v.s. time
    drive = np.zeros([4, resolution])
    for i, jj in enumerate(timeaxis):
        signal = amplitude * np.cos(jj * omega0) + dc_offset
        drive[0, i] = signal
        drive[1, i] = fitFunc_mag(signal, popt2[0], popt2[1], popt2[2], popt2[3], popt2[4])
        drive[3, i] = fitFunc_ang(signal, Ic, Cap, offset, slope)
    return drive


def get_fft_responses(drive):
    loss = np.mean(drive[1])
    # m0 = np.min(drive[0]) + (np.max(drive[0]) - np.min(drive[0])) / 2.0
    m3 = np.min(drive[3]) + (np.max(drive[3]) - np.min(drive[3])) / 2.0
    # fft_signal = np.abs(np.fft.rfft(drive[0] - m0, norm="ortho"))
    fft_mirror = np.abs(np.fft.rfft(drive[3] - m3, norm="ortho")) * loss
    return fft_mirror


fileName_ang = 'SQUID_phase.dat'
fileName_mag = 'SQUID_magnitude.dat'

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
scale2 = 1.0  # for the fitting
iguess = [Ic, Cap, offset, slope]
iguess2 = [Ic, Cap, R, offset2, scale2]
popt, pcov = curve_fit(fitFunc_ang, flux, data[2], p0=iguess)
data_mag[2] = data_mag[2] / np.max(data_mag[2])
popt2, pcov2 = curve_fit(fitFunc_mag, flux, data_mag[2], p0=iguess2, maxfev=500000, xtol=1e-36)
Ic, Cap, offset, slope = popt
print(popt2)
print(pcov2)
resolution = 2**16
pumpfreq = 8.9e9
omega0 = 2.0 * np.pi * pumpfreq
timeaxis = np.linspace(0.0, 200e-9, resolution)
freqaxis = np.fft.rfftfreq(timeaxis.shape[-1], (timeaxis[1] - timeaxis[0]))
# gap_freq = 88e9  # kept FFT nyquist limit below this (limit the FFT resolution)
pump_idx = np.argmin(abs(freqaxis - pumpfreq))
scale = 0.01
fluxpoints = 1
powerpoints = 101
# dc_offsets = np.linspace(-0.65, 0.75, fluxpoints)
dc_offsets = np.linspace(-0.45, -0.45, fluxpoints)
amplitudes = np.linspace(0.001, 0.051, powerpoints)
parabolas = np.zeros([fluxpoints, powerpoints, len(freqaxis)])
parabolas_1 = np.zeros([fluxpoints, powerpoints, len(freqaxis)]) # only with pump freq
parabolas_2 = np.zeros([fluxpoints, powerpoints, len(freqaxis)]) # only at first higher harmonic
parabolas_3 = np.zeros([fluxpoints, powerpoints, len(freqaxis)]) # only at second ..
effective_drive = np.zeros([fluxpoints, powerpoints, len(freqaxis)])



for kk, dc_offset in enumerate(dc_offsets):
    print(kk)
    for jj, amplitude in enumerate(amplitudes):
        drive = get_drive(amplitude, dc_offset, omega0, timeaxis)
        fft_mirror = get_fft_responses(drive)
        amp_d = fft_mirror[pump_idx]  # FFT amp component at drive
        # Flux parabola with all Frequencies
        parabolas[kk, jj] = get_full_parabola(freqaxis, fft_mirror, scale)
        # parabolas[kk, jj] = get_re_parabola(freqaxis, fft_mirror, scale)
        # Flux parabola only at pump
        parabolas_1[kk, jj] =  get_first_parabola(freqaxis, fft_mirror, scale, fcent=8.9e9)
        parabolas_2[kk, jj] =  get_first_parabola(freqaxis, fft_mirror, scale, fcent=8.9e9*2)
        parabolas_3[kk, jj] =  get_first_parabola(freqaxis, fft_mirror, scale, fcent=8.9e9*3)
        # parabolas2[kk, jj] = parabola(amp_d * scale, pumpfreq, freqaxis)
        effective_drive[kk, jj]= fft_mirror

header = ('Units,ufo,Frequency,' + str(freqaxis[0]) + ',' + str(freqaxis[-1]) +
          ',Pump,' + str(amplitudes[0]) + ',' + str(amplitudes[-1]) +
          ',FluxPos,' + str(dc_offsets[0]) + ','+str(dc_offsets[-1])+'')
output_file1 = 'output/mtx/sim_all_parabolas.mtx'
output_file2 = 'output/mtx/effective drive.mtx'
output_file3 = 'output/mtx/ratio_1.mtx'
output_file4 = 'output/mtx/ratio_2.mtx'
output_file5 = 'output/mtx/ratio_3.mtx'
savemtx(output_file1, parabolas, header)
savemtx(output_file2, effective_drive, header)
savemtx(output_file3, (parabolas_1)/(parabolas+1e-90), header)
savemtx(output_file4, (parabolas_2)/(parabolas+1e-90), header)
savemtx(output_file5, (parabolas_3)/(parabolas+1e-90), header)

m0 = np.min(drive[0]) + (np.max(drive[0]) - np.min(drive[0])) / 2.0
m3 = np.min(drive[3]) + (np.max(drive[3]) - np.min(drive[3])) / 2.0
gp.figure(1)
gp.s([flux, data_mag[2], fitFunc_mag(flux, popt2[0], popt2[1], popt2[2], popt2[3], popt2[4])], 'output/magnitude_fit.tmp')
gp.c('plot "output/magnitude_fit.tmp" u 1:2 w l t "Magnitude data"')
gp.c('replot "output/magnitude_fit.tmp" u 1:3 w l t "fit"')
gp.c('replot "output/magnitude_fit.tmp" u 1:(($2-$3)**2) w l t "difference^2"')
gp.figure(2)
gp.s([flux, data[2], fitFunc_ang(flux, popt[0], popt[1], popt[2], popt[3])], 'output/phase_fit.tmp')
gp.c('plot "output/phase_fit.tmp" u 1:2 w l t "Phase data"')
gp.c('replot "output/phase_fit.tmp" u 1:3 w l t "fit"')
gp.c('replot "output/phase_fit.tmp" u 1:(($2-$3)**2) w l t "difference^2"')
gp.figure(4)
gp.s([timeaxis, drive[0]-m0, drive[3]-m3], 'output/drives.tmp')
gp.c('plot "output/drives.tmp" u 1:2 w l')
gp.c('replot "output/drives.tmp" u 1:3 w l')
gp.figure(7)
gp.plot([freqaxis, fft_mirror], 'output/effective signal.tmp')
gp.figure(8)
gp.plot([freqaxis, parabolas[0, jj]], 'output/parabolas.tmp')
