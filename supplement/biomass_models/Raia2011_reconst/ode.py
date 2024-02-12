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
        dydt[V.IL13] = 0 # - v[1] - v[15]
        dydt[V.Rec] = - v[1]
        dydt[V.IL13_Rec] = + v[1] - v[7]
        dydt[V.JAK2] = - v[4] + v[6] + v[16]
        dydt[V.pJAK2] = + v[4] - v[6] - v[16]
        dydt[V.p_IL13_Rec] = + v[7] - v[14]
        dydt[V.STAT5] = - v[8] + v[17]
        dydt[V.pSTAT5] = + v[8] - v[17]
        dydt[V.CD274mRNA] = + v[9]
        dydt[V.SOCS3mRNA] = + v[10]
        dydt[V.SOCS3] = + v[12]
        dydt[V.DecoyR] = - v[15]
        dydt[V.IL13_DecoyR] = + v[15]

        return dydt


def param_values():
    """Parameter values"""
    x = [1] * C.NUM
    x[C.kf1] = 0.00342
    x[C.kr1] = 0.0
    x[C.V4] = 0.15706
    x[C.K4] = 1.0e+2
    x[C.V6] = 1
    x[C.K6] = 100
    x[C.V7] = 999.63100
    x[C.K7] = 1.0
    x[C.V8] = 0.03826e+4
    x[C.K8] = 1.0e+4
    x[C.V9] = 0.00008e+3
    x[C.K9] = 1.0e+3
    x[C.n9] = 1.0
    x[C.V10] = 0.00216e+3
    x[C.K10] = 1.0e+3
    x[C.n10] = 1.0
    x[C.kf12] = 0.3
    x[C.kf14] = 0.17292
    x[C.kf15] = 0.00012
    x[C.kr15] = 0.0
    x[C.V16] = 0.00062e+2
    x[C.K16] = 1.0e+2
    x[C.V17] = 0.00034e+3
    x[C.K17] = 1.0e+3

    return x


def initial_values():
    """Values of the initial condition"""
    y0 = [0] * V.NUM
    y0[V.IL13] = 0
    y0[V.Rec] = 1.3
    y0[V.JAK2] = 2.8
    y0[V.STAT5] = 165
    y0[V.DecoyR] = 0.34
    y0[V.SHP1] = 91

    return y0
