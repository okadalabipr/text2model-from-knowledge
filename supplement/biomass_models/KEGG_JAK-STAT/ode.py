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
        dydt[V.CNTF_CTF1_CLCF1_IL27_IL6_IL11_IL13_IL31_LIF_OSM] = 0 # - v[27]
        dydt[V.Rec] = - v[27]
        dydt[V.IL13_Rec] = + v[27] - v[28]
        dydt[V.a_Rec] = + v[28] - v[30]
        dydt[V.JAK1_JAK2_CNTFR_IL31RA_IL4R_IL6R_IL6ST_IL11RA_IL13RA1_IL13RA2_LIFR_OSMR_TYK2] = - v[29] + v[31] + v[33] + v[39]
        dydt[V.a_JAK1_JAK2_CNTFR_IL31RA_IL4R_IL6R_IL6ST_IL11RA_IL13RA1_IL13RA2_LIFR_OSMR_TYK2] = + v[29] - v[31] - v[33] - v[39]
        dydt[V.STAT] = - v[32] + v[34]
        dydt[V.a_STAT] = + v[32] - v[34]
        dydt[V.SOCSmRNA] = + v[36]
        dydt[V.SOCS4_SOCS7_SOCS1_SOCS2_SOCS3_SOCS6_SOCS5] = + v[37] - v[38]

        return dydt


def param_values():
    """Parameter values"""
    x = [1] * C.NUM
    x[C.kf30] = 0.001
    x[C.V33] = 0.001
    x[C.V36] = 0.01
    x[C.K36] = 50
    x[C.kf37] = 10
    x[C.V39] = 0.001

    return x


def initial_values():
    """Values of the initial condition"""
    y0 = [0] * V.NUM
    y0[V.CNTF_CTF1_CLCF1_IL27_IL6_IL11_IL13_IL31_LIF_OSM] = 0
    y0[V.Rec] = 1.3
    y0[V.JAK1_JAK2_CNTFR_IL31RA_IL4R_IL6R_IL6ST_IL11RA_IL13RA1_IL13RA2_LIFR_OSMR_TYK2] = 2.8
    y0[V.STAT] = 165
    y0[V.PTPN6] = 91

    return y0
