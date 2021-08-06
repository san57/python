def check_options(self, chi, **kwargs):
    """Check whether necessary parameters are define.
    Fills with default values if not

    Args:
        self (Plugin): the minimizer plugin
        chi (numpy array): Chi vector

    Return:
        the updated minimizer

    """

    # Required reduction in gradient norm
    self.zreduc = getattr(self, "zreduc", 1e-15)
    self.pevbnd = getattr(self, "pevbnd", 0.01)
    self.kvadim = chi.size
    self.kverbose = getattr(self, "kverbose", 1)
    self.ldsolve = getattr(self, "ldsolve", 1)

    # Check for missing attributes
    if not hasattr(self, "maxiter"):
        raise AttributeError(
            "maxiter is missing in the definition of the "
            "minimizer. Please check your Yaml file"
        )

    else:
        self.knevecout = self.maxiter

    return self
