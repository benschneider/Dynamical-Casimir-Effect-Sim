import warnings
import numpy as np
from math import factorial
from PyGnuplot import gp
import fit_mag_phase as fit
import logging

def xderiv(dat, dx=1.0):
    '''
    This  derivative is inaccurate as the edges.
    Calculates a 3 point derivative of a 1D matrix.
    This does not require you to shift the xaxis by one half pt.
    dx = distance between points
    what happens:
    inpu  dat = 1 2 3
        a2  = 0 1 2 3 0
          m1  = 1 2 3 - 0 1 2 = 1 1 1
          m2  = 2 3 0 - 1 2 3 = 1 1 -3
          dy  = (1 1 1 + 1 1 -3)/2 = 1 1 -1
          dy  = 1 1 1
    '''
    a2 = np.zeros(len(dat)+2)
    a2[1 :-1] = dat
    m1 = dat - a2[:-2]
    m2 = a2[2 :] - dat
    dy = (m1+m2)/2.0
    dy[0] = dy[1]
    dy[-1] = dy[-2]
    return dy/dx

def load_data(fileName):
    ''' Opens data set and retuns the data as (x, y)'''
    a = open(fileName, 'r+')
    b = a.readlines()
    data = np.zeros([2, len(b[:-3])])
    for i, e in enumerate(b[:-3]):
        strings = e.split(' ')
        data[0, i] = (float(strings[0]))  # raw xaxis
        data[1, i] = (float(strings[1]))  # raw yaxis
    return data

def get_a(x, fx, nn=50):
    ''' returns derivatives in a list up to nn order
    derivative is of fx,
    x needs to be the position of the point with equal number of points on both ends.
    i.e. for derivatives at pos 0.1
    x = np.linspace(0.1 - dx, x + dx, 1+points)
    '''
    a = []
    c = fx
    print(int(len(x)/2))
    dx = x[1]-x[0]
    for i in range(nn):
        b = xderiv(c, dx)
        c = b
        a.append(b[int(len(x)/2)])
    return a

def nout(a, ac=0.02, f=4.1e9, nmax=3, jmax=3, fd=8.9e9):
    vmc = get_part3(a, ac, nmax, jmax)
    return vmc**2*(fd-f)/(4*f)

def get_part1(jj, n):
    ''' works as expected for n=>1'''
    result = 1
    for k in range(1, jj+1):
        result = 1/(4*k*(n+k))*result
        # print(k,n, 1/result)
    print('\t 1/p3: '+ str(1/result))
    return result

def get_part2(a, ac, jmax=2, n=1):
    tb = 0
    for jj in range(jmax+1):
        p2 = a[2*jj+n-1] * ac**(2*jj+n)
        print(str(jj))
        print('\t eq: ' +str(a[2*jj+n-1])+'*AC^'+str(2*jj+n))
        p3  = get_part1(jj, n)
        print('\t p3: ' +str(p3))
        tb += p2 * p3
    return tb

def get_part3(a, ac, nmax, jmax):
    vmc = 0
    for n in range(1, nmax+1):
        p = 2**(n-1)*factorial(n-1)
        vmc += 1/p*get_part2(a, ac, jmax, n)
    return -vmc

# file_phase = '1169S11_p_4p1.dat'
# file_mag = '1169S11_m_4p1.dat'
# fit assums the first column to be flux and second data
garbage = fit.get_fitvals(file_ang='1169S11_p_4p1.dat')  # runs the fit
del garbage
fit.get_ang_resp(-0.47)  #
x = np.linspace(-0.5, -0.44, 15001)
# x = np.linspace(-0.7, 0.7, 1401)
fx = fit.get_ang_resp(x)

a = get_a(x, fx, nn=50)  # get derivatives at -0.47
logging.basicConfig(level=logging.INFO)

# fig1 = gp()
# fig1.s([*a], 'ttt.txt')
