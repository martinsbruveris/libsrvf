from scipy.interpolate import make_lsq_spline, BSpline
from interface import implements
from libsrvf.curve import Curve
import numpy as np

class BsplineLsqCurve(implements(Curve)):

    def __init__(self, t, x, degree=3, knots=None):
        """ Creates a B-spline representing a least-squares spline approximation
        to the given sample points.

        Args:
            t (array): Abscissas
            x (ndarray): The sample points. x[i] contains the coordinates of the ith
                         sample point in R^n.
            knots (array): Must satisfy the Schoenberg-Whitney conditions, i.e., there must be a subset of data points knots[j] such that knots[j] < t[j] < knots[j+k+1], for j=0, 1,...,n-k-2.
        """
        if not knots:
            knots = self.make_default_knots(t, degree)

        self.domain = (t[0], t[-1])
        self.dim = x.shape[1]
        self.f = [None] * self.dim
        self.df = [None] * self.dim
        for i in range(0, self.dim):
            xi = [xj[i] for xj in x]
            self.f[i] = make_lsq_spline(t, xi, knots, degree)
            self.df[i] = self.f[i].derivative(1)

    def make_default_knots(self, t, k):
        return np.r_[(t[0],)*k, np.linspace(t[0], t[-1], len(t)/k), (t[-1],)*k]

    def get_dimension(self):
        return self.dim

    def get_domain(self):
        return self.domain

    def evaluate(self, t):
        return [[self.f[i](tj) for i in range(self.dim)] for tj in t]

    def __call__(self, t):
        return [[self.f[i](tj) for i in range(self.dim)] for tj in t]

    def evaluate_derivative(self, t):
        return [[self.df[i](tj) for i in range(self.dim)] for tj in t]
