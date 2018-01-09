from scipy.interpolate import make_lsq_spline, BSpline
from interface import implements
from libsrvf.curve import Curve
import numpy as np

class PiecewiseLinearCurve(implements(Curve)):

    def __init__(self, parameters, samples):
        """
        Args:
            parameters (array): Parameter values corresponding to the sample points
            samples (ndarray): Sample points
        """
        self.parameters = parameters
        self.samples = samples

    def evaluate(self, t):
        pass

    def get_dimension(self):
        return self.samples.shape[1]

    def get_domain(self):
        return (self.parameters[0], self.parameters[-1])

    def evaluate(self, t):
        pass

    def evaluate_derivative(self, t):
        pass

    def find_subintervals(self, t):
        self.check_interpolation_parameters(t)
        idx = 0
        result = []
        for ti in t:
            while idx + 2 < len(self.parameters) and ti >= self.parameters[idx + 1]:
                idx += 1
            result.append(idx)
        return result

    def check_interpolation_parameters(self, t):
        if not np.all(np.diff(t) >= 0):
            raise Exception("t must be non-decreasing")
        elif t[0] < self.parameters[0] or t[-1] > self.parameters[-1]:
            raise Exception("t out of range")
