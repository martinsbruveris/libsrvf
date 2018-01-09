import numpy as np
from libsrvf.curve_impl import bspline
from libsrvf import fileio_mat
import matplotlib.pyplot as plt

def test_foo():
    path = '/home/yevgeni/Documents/research/datasets/mpeg/categorized/apple/apple-13.mat'
    mf = fileio_mat.Matfile(path)
    data = mf.get_ndarray().transpose()
    (nsamps, dim) = data.shape
    abscissas = np.linspace(0, 1, nsamps)
    f = bspline.BsplineLsqCurve(abscissas, data)
    t = np.linspace(0, 1, nsamps * 2)

#    plot2d(f(t), 'b-')
#    plt.show()

#    plot2d(f.evaluate_derivative(t), 'b-')
#    plt.show()

def plot2d(samps, color):
    sampsx = [s[0] for s in samps]
    sampsy = [s[1] for s in samps]
    plt.plot(sampsx, sampsy, color, ms=5)
