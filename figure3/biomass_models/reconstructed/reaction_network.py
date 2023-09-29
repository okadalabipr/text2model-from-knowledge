from typing import Dict, List

from .name2idx import C, V


class ReactionNetwork(object):
    """
    Reaction indices grouped according to biological processes.
    This is used for sensitivity analysis (``target``='reaction').
    """

    def __init__(self) -> None:
        super(ReactionNetwork, self).__init__()
        self.reactions: Dict[str, List[int]] = {}

    @staticmethod
    def flux(t, y, x):
        """
        Flux vector.
        """

        v = {}
        v[1] = x[C.kf1] * y[V.EGF] * y[V.EGFR] - x[C.kr1] * y[V.Ra]
        v[2] = x[C.kf2] * y[V.Ra] * y[V.Ra] - x[C.kr2] * y[V.R2]
        v[3] = x[C.kf3] * y[V.R2]
        v[4] = x[C.V4] * y[V.RP] / (x[C.K4] + y[V.RP])
        v[5] = x[C.kf5] * y[V.PLCg] * y[V.RP] - x[C.kr5] * y[V.RPL]
        v[6] = x[C.V6] * y[V.RP] * y[V.PLCg] / (x[C.K6] + y[V.PLCg])
        v[7] = x[C.kf7] * y[V.RPLP] - x[C.kr7] * y[V.RP] * y[V.PLCgP]
        v[8] = x[C.V8] * y[V.PLCgP] / (x[C.K8] + y[V.PLCgP])
        v[9] = x[C.kf9] * y[V.Grb2] * y[V.RP] - x[C.kr9] * y[V.RG]
        v[10] = x[C.kf10] * y[V.SOS] * y[V.RG] - x[C.kr10] * y[V.RGS]
        v[11] = x[C.kf11] * y[V.RGS] - x[C.kr11] * y[V.RP] * y[V.GS]
        v[13] = x[C.kf13] * y[V.Shc] * y[V.RP] - x[C.kr13] * y[V.RSH]
        v[14] = x[C.V14] * y[V.RP] * y[V.Shc] / (x[C.K14] + y[V.Shc])
        v[15] = x[C.kf15] * y[V.RShP] - x[C.kr15] * y[V.RP] * y[V.ShP]
        v[16] = x[C.V16] * y[V.ShP] / (x[C.K16] + y[V.ShP])
        v[20] = x[C.kf20] * y[V.RShGS] - x[C.kr20] * y[V.RP] * y[V.ShGS]
        v[21] = x[C.kf21] * y[V.Grb2] * y[V.ShP] - x[C.kr21] * y[V.ShG]
        v[22] = x[C.kf22] * y[V.SOS] * y[V.ShG] - x[C.kr22] * y[V.ShGS]
        v[23] = x[C.kf23] * y[V.ShGS] - x[C.kr23] * y[V.ShP] * y[V.GS]
        v[24] = x[C.kf24] * y[V.RShP] * y[V.GS] - x[C.kr24] * y[V.RShGS]
        v[25] = x[C.kf25] * y[V.PLCgP] - x[C.kr25] * y[V.PLCgP_I]
        v[28] = x[C.kf28] * y[V.GS] - x[C.kr28] * y[V.Grb2] * y[V.SOS]

        return v
