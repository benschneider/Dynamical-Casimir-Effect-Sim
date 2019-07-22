import numpy as np
from struct import pack, unpack
from codecs import open

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


def loadmtx(filename, dims=False):
    '''
    Loads an mtx file (binary compressed file)
    (first two lines of the MTX contain information of the data shape and
    what units, limits are present)
    i.e.:
    mtx, header = loadmtx('filename.mtx')
    mtx     :   will contain a 3d numpy array of the data
    header  :   will contain information on the labels and limits
    '''
    with open(filename, 'rb') as f:
        line = f.readline()
        header = line[:-1].decode('ascii')
        line = f.readline()
        a = line[:-1].decode('ascii')
        a = a.split(' ')
        s = np.array([int(strval) for strval in a])
        raw = f.read()  # reads everything else
        f.close()
    if s[3] == 4:
        data = unpack('f' * (s[2] * s[1] * s[0]), raw)  # uses float
        M = np.reshape(data, (s[2], s[1], s[0]), order="F")
    else:
        data = unpack('d' * (s[2] * s[1] * s[0]), raw)  # uses double
        M = np.reshape(data, (s[2], s[1], s[0]), order="F")
    return M, header

if __name__ == "__main__":
    x = np.linspace(-0.75,0.7011,101)
    y = np.linspace(-10, 10, 101)
    z = np.zeros([1, 101, 101])
    for i, xval in enumerate(x):
        for j, yval in enumerate(y):
            z[0, i, j] = np.float(np.sin(2*np.pi*xval)*np.cos(np.pi*2*yval))
    savemtx('test.mtx', z)

    zmat, h = loadmtx('test.mtx')
    print(zmat[0,:,0])

