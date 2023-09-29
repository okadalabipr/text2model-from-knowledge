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
            'Total_phosphorylated_Shc',
            'Total_Grb2_bound_to_Shc',
            'Total_Grb2_bound_to_EGFR',
            'Total_phosphorylated_PLCg',
            'Total_phosphorylated_EGFR',
        ]
        self.t: range = range(0, 120+1)
        self.conditions: list = [
            'EGF20nM',
            'EGF2nM',
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

        for i, condition in enumerate(self.conditions):
            if condition == 'EGF20nM':
                y0[V.EGF] = 680
            elif condition == 'EGF2nM':
                y0[V.EGF] = 68


            sol = solve_ode(self.diffeq, y0, self.t, tuple(x))

            if sol is None:
                return False
            else:
                self.simulations[self.obs_names.index('Total_phosphorylated_Shc'), i] = (
                    (sol.y[V.RShP] + sol.y[V.RShG] + sol.y[V.RShGS] + sol.y[V.ShP] + sol.y[V.ShG] + sol.y[V.ShGS]) / (sol.y[V.RShP] + sol.y[V.RShG] + sol.y[V.RShGS] + sol.y[V.ShP] + sol.y[V.ShG] + sol.y[V.ShGS] + sol.y[V.Shc] + sol.y[V.RSh])
                )
                self.simulations[self.obs_names.index('Total_Grb2_bound_to_Shc'), i] = (
                    (sol.y[V.RShG] + sol.y[V.ShG] + sol.y[V.RShGS] + sol.y[V.ShGS]) / (sol.y[V.RShG] + sol.y[V.ShG] + sol.y[V.RShGS] + sol.y[V.ShGS] + sol.y[V.Grb2] + sol.y[V.RG] + sol.y[V.RGS] + sol.y[V.GS])
                )
                self.simulations[self.obs_names.index('Total_Grb2_bound_to_EGFR'), i] = (
                    (sol.y[V.RG] + sol.y[V.RGS] + sol.y[V.RShG] + sol.y[V.RShGS]) / (sol.y[V.RG] + sol.y[V.RGS] + sol.y[V.RShG] + sol.y[V.RShGS] + sol.y[V.Grb2] + sol.y[V.GS] + sol.y[V.ShG] + sol.y[V.ShGS])
                )
                self.simulations[self.obs_names.index('Total_phosphorylated_PLCg'), i] = (
                    (sol.y[V.RPLP] + sol.y[V.PLCgP]) / (sol.y[V.RPLP] + sol.y[V.PLCgP] + sol.y[V.PLCg] + sol.y[V.RPL] + sol.y[V.PLCgP_I])
                )
                self.simulations[self.obs_names.index('Total_phosphorylated_EGFR'), i] = (
                    2 * (sol.y[V.RP] + sol.y[V.RPL] + sol.y[V.RPLP] + sol.y[V.RG] + sol.y[V.RGS] + sol.y[V.RSh] + sol.y[V.RShP] + sol.y[V.RShG] + sol.y[V.RShGS]) / (sol.y[V.EGFR] + sol.y[V.Ra] + 2 * (sol.y[V.R2] + sol.y[V.RP] + sol.y[V.RPL] + sol.y[V.RPLP] + sol.y[V.RG] + sol.y[V.RGS] + sol.y[V.RSh] + sol.y[V.RShP] + sol.y[V.RShG] + sol.y[V.RShGS]))
                )

    def set_data(self) -> None:
        self.experiments[self.obs_names.index("Total_phosphorylated_PLCg")] = {
            "EGF20nM": [0.0, 0.095, 0.055, 0.042, 0.038, 0.025],
            "EGF2nM": [0.0, 0.044, 0.048, 0.024, 0.013, 0.022]
        }
        self.experiments[self.obs_names.index("Total_phosphorylated_Shc")] = {
            "EGF20nM": [0.0, 0.288, 0.250, 0.261, 0.287, 0.283]
        }
        self.experiments[self.obs_names.index("Total_Grb2_bound_to_Shc")] = {
            "EGF20nM": [0.045, 0.350, 0.332, 0.302, 0.299, 0.291]
        }
        self.experiments[self.obs_names.index("Total_Grb2_bound_to_EGFR")] = {
            "EGF20nM": [0.0, 0.26, 0.24, 0.08, 0.055, 0.04]
        }
        self.experiments[self.obs_names.index("Total_phosphorylated_EGFR")] = {
            "EGF20nM": [0.0, 0.70, 0.48, 0.24, 0.22, 0.12],
            "EGF2nM": [0.0, 0.31, 0.24, 0.15, 0.13, 0.11]
        }

    def get_timepoint(self, obs_name: str) -> List[int]:
        if obs_name in self.obs_names:
            return [0, 15, 30, 45, 60, 120]
        # assert
