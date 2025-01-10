import numpy as np
from emcee.moves.red_blue import RedBlueMove

__all__ = ["StretchMoveInteger"]

class StretchMoveInteger(RedBlueMove):
    """
    A modified version of the Goodman & Weare (2010) "stretch move" where
    selected parameters are rounded to integers during initialization.
    """

    def __init__(self, a=2.0, integer_params=[], **kwargs):
        """
        Initialize the StretchMoveInteger class with integer rounding for specific parameters.

        :param a: (optional) The stretch scale parameter. Default: 2.0
        :param integer_params: List of parameter indices to round to integers. Default is None.
        :param kwargs: Additional parameters for the move.
        """
        self.a = a
        self.integer_params = integer_params  # List of parameter indices to round
        self.kwargs = kwargs
        super(StretchMoveInteger, self).__init__(**kwargs)

    def get_proposal(self, s, c, random):
        """
        Generate a proposal for the MCMC move, applying integer rounding for selected parameters.

        :param s: Current position (steps) in the chain.
        :param c: The chain of positions (states).
        :param random: Random number generator instance.

        :return: Proposed new positions and the factors for the proposal.
        """
        c = np.concatenate(c, axis=0)
        Ns, Nc = len(s), len(c)
        ndim = s.shape[1]

        # Stretch factor for the move
        zz = ((self.a - 1.0) * random.rand(Ns) + 1) ** 2.0 / self.a
        factors = (ndim - 1.0) * np.log(zz)

        # Select random indices for the move
        rint = random.randint(Nc, size=(Ns,))
        
        # Apply the stretch move, but round specific parameters to integers
        proposed_positions = c[rint] - (c[rint] - s) * zz[:, None]

        # Apply integer rounding for specific parameters
        for param_idx in self.integer_params:
            proposed_positions[:, param_idx] = np.round(proposed_positions[:, param_idx]/5)*5

        return proposed_positions, factors
