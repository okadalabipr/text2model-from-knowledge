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
        v[27] = x[C.kf27] * y[V.CNTF_CTF1_CLCF1_IL27_IL6_IL11_IL13_IL31_LIF_OSM] * y[V.Rec] * 2.2650 - x[C.kr27] * y[V.IL13_Rec]
        v[28] = x[C.kf28] * y[V.IL13_Rec] - x[C.kr28] * y[V.a_Rec]
        v[29] = x[C.V29] * y[V.a_Rec] * y[V.JAK1_JAK2_CNTFR_IL31RA_IL4R_IL6R_IL6ST_IL11RA_IL13RA1_IL13RA2_LIFR_OSMR_TYK2] / (x[C.K29] + y[V.JAK1_JAK2_CNTFR_IL31RA_IL4R_IL6R_IL6ST_IL11RA_IL13RA1_IL13RA2_LIFR_OSMR_TYK2])
        v[30] = x[C.kf30] * y[V.a_Rec]
        v[31] = x[C.V31] * y[V.a_JAK1_JAK2_CNTFR_IL31RA_IL4R_IL6R_IL6ST_IL11RA_IL13RA1_IL13RA2_LIFR_OSMR_TYK2] / (x[C.K31] + y[V.a_JAK1_JAK2_CNTFR_IL31RA_IL4R_IL6R_IL6ST_IL11RA_IL13RA1_IL13RA2_LIFR_OSMR_TYK2])
        v[32] = x[C.V32] * y[V.a_JAK1_JAK2_CNTFR_IL31RA_IL4R_IL6R_IL6ST_IL11RA_IL13RA1_IL13RA2_LIFR_OSMR_TYK2] * y[V.STAT] / (x[C.K32] + y[V.STAT])
        v[33] = x[C.V33] * y[V.PTPN6] * y[V.a_JAK1_JAK2_CNTFR_IL31RA_IL4R_IL6R_IL6ST_IL11RA_IL13RA1_IL13RA2_LIFR_OSMR_TYK2] / (x[C.K33] + y[V.a_JAK1_JAK2_CNTFR_IL31RA_IL4R_IL6R_IL6ST_IL11RA_IL13RA1_IL13RA2_LIFR_OSMR_TYK2])
        v[34] = x[C.V34] * y[V.a_STAT] / (x[C.K34] + y[V.a_STAT])
        v[36] = x[C.V36] * y[V.a_STAT] ** x[C.n36] / (x[C.K36] ** x[C.n36] + y[V.a_STAT] ** x[C.n36])
        v[37] = x[C.kf37] * y[V.SOCSmRNA]
        v[38] = x[C.kf38] * y[V.SOCS4_SOCS7_SOCS1_SOCS2_SOCS3_SOCS6_SOCS5]
        v[39] = x[C.V39] * y[V.SOCS4_SOCS7_SOCS1_SOCS2_SOCS3_SOCS6_SOCS5] * y[V.a_JAK1_JAK2_CNTFR_IL31RA_IL4R_IL6R_IL6ST_IL11RA_IL13RA1_IL13RA2_LIFR_OSMR_TYK2] / (x[C.K39] + y[V.a_JAK1_JAK2_CNTFR_IL31RA_IL4R_IL6R_IL6ST_IL11RA_IL13RA1_IL13RA2_LIFR_OSMR_TYK2])

        return v
