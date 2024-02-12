from typing import List

import numpy as np

from biomass.dynamics.solver import *

from .name2idx import C, V
from .ode import DifferentialEquation


class Observable(DifferentialEquation):
    """
    Correlating model simulations and experimental measurements.

    Attributes
    ----------
    obs_names : list of strings
        Names of model observables.

    t : range
        Simulation time span.

    conditions : list of strings
        Experimental conditions.

    simulations : numpy.ndarray
        The numpy array to store simulation results.

    normalization : nested dict
        * 'timepoint' : Optional[int]
            The time point at which simulated values are normalized.
            If :obj:`None`, the maximum value will be used for normalization.

        * 'condition' : list of strings
            The experimental conditions to use for normalization.
            If empty, all conditions defined in ``sim.conditions`` will be used.

    experiments : list of dict
        Time series data.

    error_bars : list of dict
        Error bars to show in figures.

    """

    def __init__(self):
        super(Observable, self).__init__(perturbation={})
        self.obs_names: list = [
            'IL13_stimulation',
            'RacSurf',
            'IL13_cell',
            'pIL4Ra',
            'pJAK2',
            'SOCS3mRNA',
            'SOCS3',
            'pSTAT5',
        ]
        self.t: range = range(0, 120+1)
        self.conditions: list = [
            'IL13_4',
            'IL13_20',
        ]
        self.simulations: np.ndarray = np.empty(
            (len(self.obs_names), len(self.conditions), len(self.t))
        )
        self.normalization: dict = {}
        self.experiments: list = [None] * len(self.obs_names)
        self.error_bars: list = [None] * len(self.obs_names)

    def simulate(self, x: list, y0: list, _perturbation: dict = {}):
        if _perturbation:
            self.perturbation = _perturbation
        # unperturbed steady state
        y0[V.CNTF_CTF1_CLCF1_IL27_IL6_IL11_IL13_IL31_LIF_OSM] = 0
        y0 = get_steady_state(self.diffeq, y0, tuple(x))
        if not y0:
            return False
        for i, condition in enumerate(self.conditions):
            if condition == 'no_stim':
                y0[V.CNTF_CTF1_CLCF1_IL27_IL6_IL11_IL13_IL31_LIF_OSM] = 0.0
            elif condition == 'IL13_4':
                y0[V.CNTF_CTF1_CLCF1_IL27_IL6_IL11_IL13_IL31_LIF_OSM] = 4.0
            elif condition == 'IL13_20':
                y0[V.CNTF_CTF1_CLCF1_IL27_IL6_IL11_IL13_IL31_LIF_OSM] = 20.0
            elif condition == 'IL13_40':
                y0[V.CNTF_CTF1_CLCF1_IL27_IL6_IL11_IL13_IL31_LIF_OSM] = 40.0
            elif condition == 'IL13_80':
                y0[V.CNTF_CTF1_CLCF1_IL27_IL6_IL11_IL13_IL31_LIF_OSM] = 80.0


            sol = solve_ode(self.diffeq, y0, self.t, tuple(x))

            if sol is None:
                return False
            else:
                self.simulations[self.obs_names.index('IL13_stimulation'), i] = (
                    sol.y[V.CNTF_CTF1_CLCF1_IL27_IL6_IL11_IL13_IL31_LIF_OSM]
                )
                self.simulations[self.obs_names.index('RacSurf'), i] = (
                    sol.y[V.Rec] + sol.y[V.IL13_Rec] + sol.y[V.a_Rec]
                )
                self.simulations[self.obs_names.index('IL13_cell'), i] = (
                    5.56750 * (sol.y[V.IL13_Rec] + sol.y[V.a_Rec])
                )
                self.simulations[self.obs_names.index('pIL4Ra'), i] = (
                    1.88700 * sol.y[V.a_Rec]
                )
                self.simulations[self.obs_names.index('pJAK2'), i] = (
                    1.39040 * sol.y[V.a_JAK1_JAK2_CNTFR_IL31RA_IL4R_IL6R_IL6ST_IL11RA_IL13RA1_IL13RA2_LIFR_OSMR_TYK2]
                )
                self.simulations[self.obs_names.index('SOCS3mRNA'), i] = (
                    17.66990 * sol.y[V.SOCSmRNA]
                )
                self.simulations[self.obs_names.index('SOCS3'), i] = (
                    sol.y[V.SOCS4_SOCS7_SOCS1_SOCS2_SOCS3_SOCS6_SOCS5]
                )
                self.simulations[self.obs_names.index('pSTAT5'), i] = (
                    sol.y[V.a_STAT]
                )

    def set_data(self) -> None:
        self.experiments[self.obs_names.index("RacSurf")] = {
            "IL13_4": [1.41, 1.37, 1.33, 1.21, 1.26, 1.26, 1.23],
        }
        self.experiments[self.obs_names.index("IL13_cell")] = {
            "IL13_4": [1.71, 2.08, 2.37, 2.37, 2.88, 3.03, 3.10, 3.24, 3.61, 3.79, 3.17],
            "IL13_20": [4.59, 4.74, 5.61, 5.54, 6.64, 6.96, 6.49, 5.69, 6.82, 6.60, 5.69],
        }
        self.experiments[self.obs_names.index("pIL4Ra")] = {
            "IL13_4": [0.34, 1.08, 0.54, 0.30, 0.84, 1.17, 1.35, 1.45, 1.12, 1.55, 1.15, 0.99, 1.25, 1.27, 1.19, 1.13],
            "IL13_20": [0.18, 0.45, 1.63, 1.92, 1.40, 2.23, 2.25, 1.54, 1.66, 2.28, 1.36]
        }
        self.experiments[self.obs_names.index("pJAK2")] = {
            "IL13_4": [0.18, 0.54, 0.64, 0.64, 0.25, 0.15, 0.83, 1.46, 1.28, 1.28, 1.43, 0.84, 1.37, 1.76, 1.44, 1.04, 1.16, 1.37, 1.22, 1.49, 0.36, 0.67],
            "IL13_20": [0.22, 0.70, 1.28, 1.68, 1.39, 2.02, 2.11, 1.65, 2.06, 1.76, 1.82, 1.10, 1.48, 1.58, 0.92, 1.21]
        }
        self.experiments[self.obs_names.index("SOCS3mRNA")] = {
            "no_stim": [0.0, 0.0, 14.49, 5.80, 8.70, 2.90, 8.70],
            "IL13_20": [2.90, 11.59, 402.90, 281.60, 637.68, 214.49, 362.32]
        }
        # self.experiments[self.obs_names.index("CD274mRNA")] = {
        #     "IL13_20": [0.99, 0.79, 1.67, 3.30, 2.77]
        # }
        self.experiments[self.obs_names.index("SOCS3")] = {
            "IL13_20": [36.42, 23.96, 78.59, 149.52, 260.70, 273.16, 254.95, 204.15],
        }
        self.experiments[self.obs_names.index("pSTAT5")] = {
            "IL13_4": [1.61, 1.61, 7.25, 30.59, 34.62, 53.94, 48.3, 57.16, 68.43, 57.16, 83.72, 80.5, 66.01, 57.96, 52.33],
            "IL13_20": [1.61, 10.47, 32.2, 50.72, 69.23, 107.07, 111.90, 106.26, 103.85, 125.58, 115.12, 127.19, 123.17, 89.36, 66.01, 107.07],
        }

    def get_timepoint(self, obs_name: str) -> List[int]:
        if obs_name in self.obs_names:
            if obs_name == "RacSurf":
                return {
                    "IL13_4": [0, 10, 20, 30, 60, 90, 120],
                    }
            elif obs_name == "IL13_cell":
                return {
                    "IL13_4": [5, 10, 15, 20, 30, 40, 50, 60, 80, 100, 120],
                    "IL13_20": [5, 10, 15, 20, 30, 40, 50, 60, 80, 100, 120],
                }
            elif obs_name == "pIL4Ra":
                return {
                    "IL13_4": [0, 4, 5, 7, 10, 15, 17.5, 20, 22.5, 25, 30, 35, 40, 45, 50, 60],
                    "IL13_20": [0, 4, 5, 7, 10, 15, 20, 30, 40, 45, 60]
                }
            elif obs_name == "pJAK2":
                return {
                    "IL13_4": [0, 2.5, 4, 5, 7, 7.5, 10, 12.5, 15, 17.5, 20, 25, 30, 35, 40, 45, 50, 60, 75, 90, 105, 120],
                    "IL13_20": [0, 4, 5, 7, 10, 15, 20, 25, 30, 40, 45, 50, 60, 80, 100, 120]
                }
            elif obs_name == "SOCS3mRNA":
                return {
                    "IL13_20": [0, 20, 40, 60, 80, 100, 120]
                }
            # elif obs_name == "CD274mRNA":
            #     return {
            #         "IL13_20": [0, 30, 60, 90, 120]
            #     }
            elif obs_name == "SOCS3":
                return {
                    "IL13_20": [0, 20, 40, 60, 80, 90, 100, 120],
                }
            elif obs_name == "pSTAT5":
                return {
                    "IL13_4": [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 100],
                    "IL13_20": [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100]
                }
        assert False
