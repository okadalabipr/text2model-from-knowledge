import numpy as np

from biomass.estimation import convert_scale, initialize_search_param

from .name2idx import C, V
from .ode import initial_values, param_values


class SearchParam(object):
    """Specify model parameters and/or initial values to optimize."""

    def __init__(self):
        # parameters
        self.idx_params = [
            C.kf27,
            C.kr27,
            C.kf28,
            C.kr28,
            C.V29,
            C.K29,
            C.kf30,
            C.V31,
            C.K31,
            C.V32,
            C.K32,
            C.V33,
            C.K33,
            C.V34,
            C.K34,
            C.V36,
            C.K36,
            C.n36,
            C.kf37,
            C.kf38,
            C.V39,
            C.K39,
        ]

        # initial values
        self.idx_initials = []

    def get_region(self):
        x = param_values()
        y0 = initial_values()

        search_param = initialize_search_param(
            parameters=C.NAMES,
            species=V.NAMES,
            param_values=x,
            initial_values=y0,
            estimated_params=self.idx_params,
            estimated_initials=self.idx_initials,
        )

        search_rgn = np.zeros((2, len(x) + len(y0)))
        # Default: 0.1 ~ 10
        for i, j in enumerate(self.idx_params):
            search_rgn[0, j] = search_param[i] * 0.1  # lower bound
            search_rgn[1, j] = search_param[i] * 10.0  # upper bound
        # Default: 0.5 ~ 2
        for i, j in enumerate(self.idx_initials):
            search_rgn[0, j + len(x)] = search_param[i + len(self.idx_params)] * 0.5  # lower bound
            search_rgn[1, j + len(x)] = search_param[i + len(self.idx_params)] * 2.0  # upper bound

        # search_rgn[:,C.parameter] = [lower_bound, upper_bound]
        # search_rgn[:,V.specie+len(x)] = [lower_bound, upper_bound]

        search_rgn = convert_scale(
            region=search_rgn,
            parameters=C.NAMES,
            species=V.NAMES,
            estimated_params=self.idx_params,
            estimated_initials=self.idx_initials,
        )

        return search_rgn

    def update(self, indiv):
        x = param_values()
        y0 = initial_values()

        for i, j in enumerate(self.idx_params):
            x[j] = indiv[i]
        for i, j in enumerate(self.idx_initials):
            y0[j] = indiv[i + len(self.idx_params)]

        # parameter constraints
        

        return x, y0
