from interface import Interface

class Curve(Interface):
    """Represents an absolutely continuous function from an interval to R^n, n > 0.

    This is the main type that API functions take as input.  Clients of the library
    create classes that implement this interface.

    A key point is that computation of derivatives is left to the client.
    """

    def get_dimension(self):
        """Returns the dimension of the codomain.

        Returns:
            A positive integer n, representing the number of component functions.
        """
        pass

    def get_domain(self):
        """Returns the domain of the function.

        Returns:
            A 2-tuple, representing the lower and upper bounds of the domain interval.
        """
        pass

    def evaluate(self, t):
        """Evaluate the function at a point in its domain.

        Args:
            t (float): The point at which to evaluate the function.  Must lie within
                the function's domain.

        Returns:
            An n-tuple of floats, where n = self.get_dimension(), representing the
            value of the function at t.
        """
        pass

    def evaluate_derivative(self, t):
        """Evaluate the first derivative of the function at a point in its domain.

        Args:
            t (float): The point at which to evaluate the derivative.  Must lie within
                the function's domain.

        Returns:
            An n-tuple of floats, where n = self.get_dimension(), representing the value
            of the first derivative at t.
        """
        pass
