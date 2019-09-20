from scipy.optimize import curve_fit, minimize
from scipy import constants as cons
import numpy as np
from numpy import pi, cos, sin
import PyGnuplot as gp
from struct import pack, unpack
import pyfftw
from time import time
pyfftw.interfaces.cache.enable() # enable cash -> speedier execution

# import fit_mag_phase as fit_data

def get_drive(ac_offset, dc_offset, f, timeaxis, fit_data, len_drive):
    # 0: Input amplitude v.s. time
    # 1: Corresponding phase v.s. time
    drive = np.zeros([2, len_drive])  # 7us
    input_drives = ac_offset * cos(timeaxis * 2*pi*f) + dc_offset  # 0.5ms
    phase_res = fit_data.get_ang_resp(input_drives)  # 19ms
    drive[1, :] = phase_res
    drive[0, :] = input_drives
    return drive

def get_fft_response(drive, len_drive):
    # fft_mirror = np.abs(np.fft.rfft(drive[1], norm="ortho")) # use timestep not 1/sqrt(n)
    if len_drive < 10000:
        fft_mirror = np.fft.rfft(drive[1])
    elif len_drive <60000:
        fft_mirror = pyfftw.interfaces.numpy_fft.rfft(drive[1], threads=2)
    else:
        fft_mirror = pyfftw.interfaces.numpy_fft.rfft(drive[1], threads=4)
    return 2*np.abs(fft_mirror)/len_drive

def extract_hoa_values(fft_mirror, faxis, drive_freq, harmonics = 6):
    # no need to store everything just the main harmonics
    # harmonics = 10 , number of harmonics to go after
    width = 4 # integration range / range to find peak value for
    hao_ax = np.zeros(harmonics)
    for i in range(harmonics):
        pos = np.argmin(np.abs(faxis-drive_freq*(i+1)))
        for j in range(width+1):
            intval = abs(int(pos+j-width/2))
            hao_ax[i] += abs(fft_mirror[intval])
        # hao_ax[i] = fft_mirror[pos]
    return hao_ax

def make_response(ac_flux, dc_flux, drive_freq, timeaxis, fit_data, harmonics=6, plot=False):
    time_step = (timeaxis[1]-timeaxis[0])
    faxis = np.fft.rfftfreq(len(timeaxis), time_step)
    hoa_data = np.zeros([len(dc_flux), len(ac_flux), harmonics])
    len_drive = len(timeaxis)
    for ii, dc_offset in enumerate(dc_flux):
        print(ii)
        for jj, ac_offset in enumerate(ac_flux):
            drive = get_drive(ac_offset, dc_offset, drive_freq, timeaxis, fit_data, len_drive)  # 3ms
            fft_mirror = get_fft_response(drive, len_drive)  # 2ms
            hoa_ax = extract_hoa_values(fft_mirror, faxis, drive_freq, harmonics)  # 1ms
            hoa_data[ii, jj] = hoa_ax  # 1us
            # if plot:
            #     gp.figure(1)
            #     gp.plot([timeaxis,  drive[1]], 'tmp1.dat')
            #     gp.figure(2)
            #     gp.plot([faxis[:], fft_mirror[:]], 'tmp2.dat')
    return hoa_data

if __name__ == "__main__":
    #How to test feasability of functions here:
    import fit_mag_phase as fit_data
    fit_data.get_fitvals()
    x = np.linspace(-1,1,201)
    y = fit_data.get_ang_resp(x)
    # gp.plot([x,y])
    timeaxis = np.linspace(0, 5e-9, 50001)
    # faxis = np.fft.rfftfreq(len(timeaxis), (timeaxis[1]-timeaxis[0]))
    drive_freq = 8.9e9
    harmonics = 20
    dc_flux = [0.0, 1.1]
    ac_flux = [50e-3]
    matrix = make_response(ac_flux, dc_flux, drive_freq, timeaxis, fit_data, harmonics=harmonics, plot=False)
    print(matrix)

    # drive = get_drive(1, 0, drive_freq, timeaxis)
    # fft_mirror = get_fft_response(drive)
    # hoa_ax = extract_hoa_values(fft_mirror, faxis, drive_freq, harmonics)
    # print(hoa_ax)
    # gp.figure(1)
    # gp.plot([timeaxis,  drive[1]], 'tmp1.dat')
    # gp.figure(2)
    # gp.plot([faxis[10:], fft_mirror[10:]], 'tmp2.dat')


