import numpy as np

from biomass.estimation import convert_scale, initialize_search_param

from .name2idx import C, V
from .ode import initial_values, param_values


class SearchParam(object):
    """Specify model parameters and/or initial values to optimize."""

    def __init__(self):
        # parameters
        self.idx_params = [
            C.kf1,
            C.kr1,
            C.kf2,
            C.kr2,
            C.kf3,
            C.kr3,
            C.V4,
            C.K4,
            C.kf5,
            C.kr5,
            C.V6,
            C.K6,
            C.kf7,
            C.kr7,
            C.V8,
            C.K8,
            C.kf9,
            C.kr9,
            C.kf10,
            C.kr10,
            C.kf11,
            C.kr11,
            C.kf13,
            C.kr13,
            C.V14,
            C.K14,
            C.kf15,
            C.kr15,
            C.V16,
            C.K16,
            C.kf20,
            C.kr20,
            C.kf21,
            C.kr21,
            C.kf22,
            C.kr22,
            C.kf23,
            C.kr23,
            C.kf24,
            C.kr24,
            C.kf25,
            C.kr25,
            C.kf28,
            C.kr28,
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
