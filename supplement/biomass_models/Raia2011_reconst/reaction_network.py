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
        IL13_scale = 2.265
        v[1] = x[C.kf1] * y[V.IL13] * y[V.Rec] * IL13_scale - x[C.kr1] * y[V.IL13_Rec]
        v[4] = x[C.V4] * y[V.IL13_Rec] * y[V.JAK2] / (x[C.K4] + y[V.JAK2])
        v[6] = x[C.V6] * y[V.SOCS3] * y[V.pJAK2] / (x[C.K6] + y[V.pJAK2])
        v[7] = x[C.V7] * y[V.pJAK2] * y[V.IL13_Rec] / (x[C.K7] + y[V.IL13_Rec])
        v[8] = x[C.V8] * y[V.pJAK2] * y[V.STAT5] / (x[C.K8] + y[V.STAT5])
        v[9] = x[C.V9] * y[V.pSTAT5] ** x[C.n9] / (x[C.K9] ** x[C.n9] + y[V.pSTAT5] ** x[C.n9])
        v[10] = x[C.V10] * y[V.pSTAT5] ** x[C.n10] / (x[C.K10] ** x[C.n10] + y[V.pSTAT5] ** x[C.n10])
        v[12] = x[C.kf12] * y[V.SOCS3mRNA]
        v[14] = x[C.kf14] * y[V.p_IL13_Rec]
        v[15] = x[C.kf15] * y[V.IL13] * y[V.DecoyR] * IL13_scale - x[C.kr15] * y[V.IL13_DecoyR]
        v[16] = x[C.V16] * y[V.SHP1] * y[V.pJAK2] / (x[C.K16] + y[V.pJAK2])
        v[17] = x[C.V17] * y[V.SHP1] * y[V.pSTAT5] / (x[C.K17] + y[V.pSTAT5])

        return v
