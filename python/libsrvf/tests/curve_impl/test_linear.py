from libsrvf.curve_impl import linear
import numpy as np

def test_find_subintervals():
    params = [0, 0.5, 1]
    samps = np.array([[0], [0.5], [1]])
    c = linear.PiecewiseLinearCurve(params, samps)
    print c.find_subintervals([0, 0.5, 1])
    raise "foo"
