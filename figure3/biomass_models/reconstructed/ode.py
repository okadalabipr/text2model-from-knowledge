from .name2idx import C, V
from .reaction_network import ReactionNetwork


class DifferentialEquation(ReactionNetwork):
    def __init__(self, perturbation):
        super(DifferentialEquation, self).__init__()
        self.perturbation = perturbation

    def diffeq(self, t, y, *x):
        """Kinetic equations"""
        v = self.flux(t, y, x)

        if self.perturbation:
            for i, dv in self.perturbation.items():
                v[i] = v[i] * dv

        dydt = [0] * V.NUM
        dydt[V.EGF] = - v[1]
        dydt[V.EGFR] = - v[1]
        dydt[V.Ra] = + v[1] - 2 * v[2]
        dydt[V.R2] = + v[2] - v[3] + v[4]
        dydt[V.RP] = + v[3] - v[4] - v[5] + v[7] - v[9] + v[11] - v[13] + v[15] + v[20]
        dydt[V.PLCg] = - v[5] - v[6] + v[8]
        dydt[V.RPL] = + v[5]
        dydt[V.RPLP] = + v[6] - v[7]
        dydt[V.PLCgP] = + v[7] - v[8] - v[25]
        dydt[V.Grb2] = - v[9] - v[21] + v[28]
        dydt[V.RG] = + v[9] - v[10]
        dydt[V.SOS] = - v[10] - v[22] + v[28]
        dydt[V.RGS] = + v[10] - v[11]
        dydt[V.GS] = + v[11] + v[23] - v[24] - v[28]
        dydt[V.Shc] = - v[13] - v[14] + v[16]
        dydt[V.RSH] = + v[13]
        dydt[V.RShP] = + v[14] - v[15] - v[24]
        dydt[V.ShP] = + v[15] - v[16] - v[21] + v[23]
        dydt[V.RShGS] = - v[20] + v[24]
        dydt[V.ShGS] = + v[20] + v[22] - v[23]
        dydt[V.ShG] = + v[21] - v[22]
        dydt[V.PLCgP_I] = + v[25]

        return dydt


def param_values():
    """Parameter values"""
    x = [1] * C.NUM
    x[C.kf1] = 0.003
    x[C.kr1] = 0.06
    x[C.kf2] = 0.01
    x[C.kr2] = 0.1
    x[C.kf3] = 1
    x[C.kr3] = 0.01
    x[C.V4] = 450
    x[C.K4] = 50
    x[C.kf5] = 0.06
    x[C.kr5] = 0.2
    x[C.V6] = 1
    x[C.K6] = 100
    x[C.kf7] = 0.3
    x[C.kr7] = 0.006
    x[C.V8] = 1
    x[C.K8] = 100
    x[C.kf9] = 0.003
    x[C.kr9] = 0.05
    x[C.kf10] = 0.01
    x[C.kr10] = 0.06
    x[C.kf11] = 0.03
    x[C.kr11] = 4.5e-3
    x[C.kf13] = 0.09
    x[C.kr13] = 0.6
    x[C.V14] = 1
    x[C.K14] = 100
    x[C.kf15] = 0.3
    x[C.kr15] = 9e-4
    x[C.V16] = 1.7
    x[C.K16] = 340
    x[C.kf20] = 0.12
    x[C.kr20] = 2.4e-4
    x[C.kf21] = 0.003
    x[C.kr21] = 0.1
    x[C.kf22] = 0.03
    x[C.kr22] = 0.064
    x[C.kf23] = 0.1
    x[C.kr23] = 0.021
    x[C.kf24] = 0.009
    x[C.kr24] = 4.29e-2
    x[C.kf25] = 1
    x[C.kr25] = 0.03
    x[C.kf28] = 1.5e-3
    x[C.kr28] = 1e-4

    return x


def initial_values():
    """Values of the initial condition"""
    y0 = [0] * V.NUM
    y0[V.EGFR] = 100
    y0[V.PLCg] = 105
    y0[V.Grb2] = 85
    y0[V.SOS] = 34
    y0[V.Shc] = 150

    return y0
